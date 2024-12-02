# library(shiny)
# library(GGally)
# library(ggplot2)
# library(shinydashboard)

library(shinyFeedback)
library(shiny)
library(shinyjs)
library(shinydashboard)
library(DT)

# Define UI for application that draws a histogram
function(request) {
    shinyjs::useShinyjs()  # Enable shinyjs
    shinyFeedback::useShinyFeedback()  # Enable shinyFeedback
    
    
    sidebar <- dashboardSidebar(
      hr(),
      sidebarMenu(id="tabs",
                  menuItem("Analysis", tabName="Analysis", icon=icon("line-chart")),
                  menuItem("Data", tabName = "Data", icon=icon("table")),
                  menuItem("Codes",  icon = icon("file-text"),
                           menuSubItem("ui.R", tabName = "ui", icon = icon("angle-right")),
                           menuSubItem("server.R", tabName = "server", icon = icon("angle-right")),
                           menuSubItem("Functions.R", tabName = "functions", icon = icon("gears"))
                  ),
                  menuItem("Read Me", tabName = "readme", icon=icon("mortar-board"), selected=TRUE),
                  menuItem("About", tabName = "about", icon = icon("question"))
      )
    )
    
    body <- dashboardBody(
      tabItems(
        tabItem(tabName = "readme",
                fluidPage(
                  tags$iframe(src = './readme.html', 
                              width = '100%', height = '800px', 
                              frameborder = 0, scrolling = 'auto'
                  )
                )
        ),
        tabItem(tabName = "Analysis",
                fluidRow(
                  column(width = 6, 
                         tabBox(width = NULL,
                                tabPanel(
                                  h5("Simulation Parameters"),
                                  fluidRow(
                                    column(
                                      width = 12,
                                      box(
                                        width = NULL, collapsible = TRUE,
                                        title = "Cluster Info", solidHeader = TRUE,
                                        splitLayout(
                                          numericInput("G", "Number of Clusters:", value = 2, min = 1, step = 1),
                                          numericInput("d", "Number of Dimensions:", value = 3, min = 2, step = 1),
                                          numericInput("N", "Sample Size:", value = 250, min = 50, step = 10)
                                        )
                                      ))
                                  )
                                ),
                                tabPanel(
                                  h5("Cluster Parameters"),
                                  fluidRow(
                                    column(
                                      width = 12,
                                      uiOutput("cluster_parameters_ui")  # Dynamically generated UI
                                    )
                                  ),
                                  fluidRow(
                                    column(
                                      width = 12,
                                      actionButton("generate_data", "Generate Data")  # Button to trigger data generation
                                    )
                                  )
                                ),
                                tabPanel(
                                  h5("Amputation"),
                                  fluidRow(
                                    column(
                                      width = 12,
                                      box(
                                        width = NULL, collapsible = TRUE,
                                        title = "Missingness Settings", solidHeader = TRUE,
                                        checkboxInput(
                                          inputId = "ampute_data",
                                          label = "Ampute the Data?",
                                          value = FALSE
                                        ),
                                        conditionalPanel(
                                          condition = "input.ampute_data == true",
                                          splitLayout(
                                            numericInput(
                                              inputId = "percent_missing",
                                              label = "Percent of Data to Make Missing:",
                                              value = 0.10,
                                              min = 0,
                                              max = 0.80,
                                              step = 0.01
                                            ),
                                            radioButtons(
                                              inputId = "missing_mechanism",
                                              label = "Missingness Mechanism:",
                                              choices = c(
                                                "MCAR - Missing completely at random",
                                                "MNAR - Missing not at random",
                                                "MAR - Missing at random"
                                              ),
                                              selected = "MCAR - Missing completely at random",  # Set a default choice
                                              inline = FALSE  # Set to TRUE if you want horizontal alignment
                                            )
                                          ),
                                          actionButton(
                                            inputId = "ampute_button",
                                            label = "Ampute!"
                                          )
                                        )
                                      )
                                    )
                                  )
                                ),tabPanel(
                                  h5("Clustering Method"),
                                  fluidRow(
                                    column(
                                      width = 12,
                                      radioButtons(
                                        inputId = "clustering_method",
                                        label = "Choose the Clustering Method:",
                                        choices = c("MixPLN - Mixtures of Poisson log-normal",
                                                    "KMM - k missing means"),
                                        selected = "MixPLN - Mixtures of Poisson log-normal",
                                        inline = TRUE
                                      ),
                                      actionButton(
                                        inputId = "cluster_button",
                                        label = "Cluster data using selected method!"
                                      ),
                                      uiOutput("cluster_info_boxes")
                                    )
                                  )
                                ),
                                tabPanel(
                                  h5("Plotting Setting"),
                                  fluidRow(
                                    column(
                                      width = 12,
                                      uiOutput("color_selection_ui")  # Dynamically generated UI
                                    )
                                  ),
                                  fluidRow(
                                    column(
                                      width = 12,
                                      actionButton("render_plot", "Render Plot")  # Button to trigger plot rendering
                                    )
                                  )
                                )
                         )
                  ),
                  column(
                    width = 6,
                    box(
                      width = NULL,
                      plotOutput("plot1", height = "650px", brush = "plot_brush"),
                      title = "Pairs Plot of all Clusters",
                      solidHeader = TRUE,
                      status = "primary"
                    )
                  )
                )
        ),
        tabItem(
          tabName = "Data",
          box(
            width = NULL, status = "primary", solidHeader = TRUE, title = "Data Table",
            actionButton("use_example", "Load Example Dataset"),
            br(),
            DTOutput("data_table")
          )
        ),
        tabItem(tabName = "ui",
                box( width = NULL, status = "primary", solidHeader = TRUE, title="ui.R",
                     br(),br(),
                     pre(includeText("ui.R"))
                )
        ),
        tabItem(tabName = "server",
                box( width = NULL, status = "primary", solidHeader = TRUE, title="server.R",
                     br(),br(),
                     pre(includeText("server.R"))
                )
        ),
        tabItem(tabName = "about",
                fluidPage(
                  tags$iframe(src = './about.html', 
                              width = '100%', height = '800px',
                              frameborder = 0, scrolling = 'auto'
                  )
                )
                
        
        ),
        tabItem(tabName = "functions",
                box( width = NULL, status = "primary", solidHeader = TRUE, title="Functions.R",
                     br(),br(),
                     pre(includeText("Functions.R"))
                )
        )
      )
    )
    
    
    dashboardPage(
      dashboardHeader(title = "Missing Data"),
      sidebar,
      body
    )

}