# library(shiny)
# library(shinyjs)
# library(shinyFeedback)

library(shinyFeedback)
library(shiny)
library(shinyjs)
library(shinydashboard)
library(mclust)
library(GGally)
library("clusterGeneration")
library(mvtnorm)
library("mice")
library(gtools)
library(dplyr)
library(data.table)
library(DT)
library(e1071)
library(glue)
library(DangLoop)
library(kmmeans)

cluster_colors <- c("green","blue")
already_amputed <- reactiveVal(FALSE)

server <- function(input, output, session) {
  load("Example_Dataset.Rdata")  
  load("Functions.Rdata")
  
  data_reactive <- reactiveVal(dat)
  table_data <- NULL
  
  #load example data set if check box is checked
  observeEvent(c(input$use_example,table_data), {
    already_amputed(FALSE)
    load("Example_Dataset.Rdata")
    data_reactive(dat)
    print("LOADING EXAMPLE DATA")
    table_data <- data.table(data_reactive()[[2]])
    output$data_table <- DT::renderDT(
      table_data  # Display the second element of the dataset
    )
  })
 
  
  #color choices for clusters
  output$color_selection_ui <- renderUI({
    req(input$G)
    validate(need(input$G > 0, "Number of clusters must be greater than 0"))
    
    color_choices <- c("red", "green", "blue", "yellow", "purple", "orange", "cyan",
                       "pink", "brown", "tan", "grey", "black")
    
    lapply(1:input$G, function(i) {
      selectInput(
        inputId = paste0("cluster_color_", i),
        label = paste("Color for Cluster", i),
        choices = color_choices,
        selected = color_choices[i %% length(color_choices) + 1]
      )
    })
  })
  
  #obtaining parameters for cluster generation
  output$cluster_parameters_ui <- renderUI({
    req(input$G, input$G)
    validate(need(input$G > 0, "Number of clusters must be greater than 0"))
    validate(need(input$d > 0, "Number of dimensions must be greater than 0"))
    
    lapply(1:input$G, function(i) {
      tagList(
        h4(paste("Cluster", i, "Parameters")),
        lapply(1:input$d, function(j) {
          numericInput(
            inputId = paste0("cluster_", i, "_mean_dim_", j),
            label = paste("Mean for Dimension", j, "(Cluster", i, ")"),
            value = 1, min = 1, max = 10, step = 0.1
          )
        }),
        numericInput(
          inputId = paste0("cluster_", i, "_alpha"),
          label = paste("Alpha Value for Cluster", i,"(Controls covaraince structure. Larger alpha = less correlation)"),
          value = 1, min = 1, step = 0.1
        ),
        numericInput(
          inputId = paste0("cluster_", i, "_sd_lower"),
          label = paste("Lower Bound for Standard Deviation (Cluster", i, ")"),
          value = 1, min = 1, max = 6, step = 0.1
        ),
        numericInput(
          inputId = paste0("cluster_", i, "_sd_upper"),
          label = paste("Upper Bound for Standard Deviation (Cluster", i, ")"),
          value = 1.5, min = 1, max = 6, step = 0.1
        ),
        numericInput(
          inputId = paste0("cluster_", i, "_percent"),
          label = paste("Percent of Data for Cluster", i),
          value = 0.4, min = 0, step = 0.1
        ),
        hr()
      )
    })
  })
  
  
  # Reactive to extract parameters for data generation
  observeEvent(input$generate_data, {
    req(input$G, input$d)  # Ensure inputs are defined
    load("Functions.Rdata")
    validate(
      need(input$G > 0, "Number of clusters must be greater than 0"),
      need(input$d > 0, "Number of dimensions must be greater than 0")
    )
    
    parameters <- lapply(1:input$G, function(i) {
      list(
        alpha = input[[paste0("cluster_", i, "_alpha")]],
        sd_lower = input[[paste0("cluster_", i, "_sd_lower")]],
        sd_upper = input[[paste0("cluster_", i, "_sd_upper")]],
        mu = sapply(1:input$d, function(j) {
          input[[paste0("cluster_", i, "_mean_dim_", j)]]
        }),
        percent = input[[paste0("cluster_", i, "_percent")]]
      )
    })
    
    
    validate(
      need(sum(pi_g) == 1, "Cluster percentages must sum to 1")
    )
    
    #format parameters into nice format
    par <- format_params(parameters)
    
    #generate data
    dat <- generate_data(par,input$N)
    already_amputed(FALSE)
    data_reactive(dat)
    #return generated data
    print(data_reactive()[[2]])
    print("ATTEMPTING TO RENDER")
    table_data <- data.table(data_reactive()[[2]])
    output$data_table <- DT::renderDT(
      table_data  
    )
    
    output$plot1 <- renderPlot({
      plot_data <- plot_trigger()$plot_data  
      cluster_colors <- plot_trigger()$cluster_colors 
      
      plot_data[which(is.na(plot_data))] <- -666
      
      if(any(is.na(plot_data))){
        print("NA VLAUES IN PLOT DATA")
      }
      
      # Pass cluster_colors to scat.my via GGally's wrap
      ggpairs(
        plot_data,
        aes(color = Cluster),  
        lower = list(continuous = wrap(scat.my, colors = cluster_colors)),  
        upper = list(continuous = wrap("cor")), 
        diag = list(continuous = wrap(diag.my, colors = cluster_colors)
        )) +
        scale_color_manual(values = cluster_colors) + 
        theme_minimal()
    })
    
    print(data_reactive())
  })
  
  #ampute data
  observeEvent(input$ampute_button, {
    print("Ampute button pressed!")
    
    if (!already_amputed()) {  # Only ampute if not already amputed
      print("NOT AMPUTED YETT")
      req(data_reactive(), input$percent_missing, input$missing_mechanism)
      dat <- data_reactive()
      # Ampute the data
      print("ABOUT TO AMPUTRE")
      amputed_data <- generate_missing(
        dat,
        percent_miss = input$percent_missing,
        mech = input$missing_mechanism
      )
      print("EMMEMEME")
      dat <- amputed_data
      print("SETTING VALUES TO DATA_RECATIVE")
      data_reactive(dat)  
      print("AMPUTED NOW")
      already_amputed(TRUE) 
      
      output$plot1 <- renderPlot({
        plot_data <- plot_trigger()$plot_data  
        cluster_colors <- plot_trigger()$cluster_colors 
        
        plot_data[which(is.na(plot_data))] <- -666
        
        if(any(is.na(plot_data))){
          print("NA VLAUES IN PLOT DATA")
        }
        
        # Pass cluster_colors to scat.my via GGally's wrap
        ggpairs(
          plot_data,
          aes(color = Cluster),  
          lower = list(continuous = wrap(scat.my, colors = cluster_colors)),  # Customize lower panels
          upper = list(continuous = wrap("cor")), 
          diag = list(continuous = wrap(diag.my, colors = cluster_colors)
          )) +
          scale_color_manual(values = cluster_colors) + 
          theme_minimal()
      })
      
      table_data <- data.table(data_reactive()[[2]])
      output$data_table <- DT::renderDT(
        table_data
      )
      
    }else{
      print("ALREADY AMPUTED")
      shinyjs::runjs('alert("Data has already been amputed. Cannot ampute again.");')
    }

    
  })
  
  #check for when plot should be rendered
  plot_trigger <- eventReactive({c(input$render_plot,input$generate_data)}, {
    req(data_reactive(), input$G)  # Ensure data and cluster count are available
    
    data <- data_reactive()
    log_data <- log(data[[2]] + 1 / 6)  
    
    print("PASSED LOG")
    
    cluster_labels <- mclust::map(data[[4]])
    
    validate(
      need(length(cluster_labels) == nrow(data[[2]]), "Cluster labels must match the number of rows in the data")
    )
    
    plot_data <- data.frame(log_data, Cluster = factor(cluster_labels))
    
    # Get cluster colors
    cluster_colors <- sapply(1:input$G, function(i) {
      input[[paste0("cluster_color_", i)]] %||% "gray"  # Default to "gray" if a color is missing
    })
    
    validate(
      need(length(cluster_colors) >= length(unique(cluster_labels)), 
           "Not enough colors provided for the number of clusters.")
    )
    print("PASSED TRIGGER")
    list(plot_data = plot_data, cluster_colors = cluster_colors)  
  })
  
  #pairs plot for visualizing data
  output$plot1 <- renderPlot({
    plot_data <- plot_trigger()$plot_data  
    cluster_colors <- plot_trigger()$cluster_colors  
    
    plot_data[which(is.na(plot_data))] <- -666
    
    if(any(is.na(plot_data))){
      print("NA VLAUES IN PLOT DATA")
    }
    
    # Pass cluster_colors to scat.my via GGally's wrap
    ggpairs(
      plot_data,
      aes(color = Cluster),  
      lower = list(continuous = wrap(scat.my, colors = cluster_colors)),  
      upper = list(continuous = wrap("cor")),  
      diag = list(continuous = wrap(diag.my, colors = cluster_colors)
    )) +
      scale_color_manual(values = cluster_colors) + 
      theme_minimal()
  })
  
  
  cluster_outputs <- reactiveValues()
  observeEvent(input$cluster_button, {
    print("Cluster button pressed!")
    req(input$G, input$clustering_method)
    
   
    max_clusters <- input$G + 2
    
 
    cluster_results <- list()
    clustering_data <- data_reactive()
    
    print("CHECKING METHOD")
    if (grepl("MixPLN",input$clustering_method,fixed=TRUE)) {
      print("STARTING CLUSTER")
      
      
      for (cluster in 1:max_clusters) {
        cluster_results[[cluster]] <- mixPLN(clustering_data, cluster)
      }
      print("FINISH CLUSTER")
      
      
      for (cluster in 1:max_clusters) {
        cluster_results[[cluster]]$ARI <- mclust::adjustedRandIndex(
          mclust::map(cluster_results[[cluster]]$z),
          mclust::map(cluster_results[[cluster]]$true)
        )
      }
      BIC_vector <- c()
      for(cluster in 1:max_clusters){
        BIC_vector[cluster] <- cluster_results[[cluster]]$BIC
      }
      
      print("1")
      
      optimal_cluster <- which.max(BIC_vector)
      optimal_data <- cluster_results[[optimal_cluster]]
      optimal_loglik <- optimal_data$loglik
      optimal_errors <- mclust::classError(
        mclust::map(optimal_data$z),
        mclust::map(optimal_data$true)
      )
      print("2")
      optimal_mu_est <- optimal_data$mu
      optimal_BIC <- optimal_data$BIC
      optimal_ARI <- optimal_data$ARI
      
     
      
      
      
      
      output$cluster_info_boxes <- renderUI({
        req(cluster_results)  # Ensure cluster_results is available
        
        
        cluster_data <- isolate({
          lapply(1:max_clusters, function(k) {
            ARI <- cluster_results[[k]]$ARI
            BIC <- cluster_results[[k]]$BIC
            misclasss <- mclust::classError(
              mclust::map(cluster_results[[k]]$z),
              mclust::map(cluster_results[[k]]$true)
            )
            list(
              BIC = BIC,
              ARI = ARI,
              num_misclass = length(misclasss$misclassified),
              which_misclass = misclasss$misclassified,
              loglik = cluster_results[[k]]$loglik,
              mu = cluster_results[[k]]$mu
            )
          })
        })
        # Render the plots dynamically
        for (k in 1:max_clusters) {
          local_k <- k  # Ensure proper scoping
          output[[paste0("loglik_plot_", local_k)]] <- renderPlot({
            plot_loglikelihood(cluster_results[[local_k]]$loglik, local_k)
          })
        }
        # Generate the UI output
        tagList(
          lapply(1:1, function(fgj){
            data <- cluster_data[[optimal_cluster]]
            content <- glue::glue("
        <h5>Optimal Number of Clusters (by BIC):</h5>
        <p>{optimal_cluster}</p>
        <h5>Largest Bayesian Information Criterion (BIC):</h5>
        <p>{data$BIC}</p>
        <h5>Adjusted Rand Index (ARI):</h5>
        <p>{data$ARI}</p>
        <h5>Number of Misclassifications:</h5>
        <p>{data$num_misclass}</p>
        <h5>Observations Misclassified:</h5>
        <pre>{paste(data$which_misclass, collapse = ', ')}</pre>
      ")
            box(
              title = paste("Summary of Optimal Clustering Information"),
              collapsible = TRUE,
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              div(
                style = 'overflow-y: auto; max-height: 400px;',
                HTML(content)
              )
            )
        }),
          lapply(1:max_clusters, function(k) {
            data <- cluster_data[[k]]
            content <- glue::glue("
        <h5>Bayesian Information Criterion (BIC):</h5>
        <p>{data$BIC}</p>
        <h5>Adjusted Rand Index (ARI):</h5>
        <p>{data$ARI}</p>
        <h5>Number of Misclassifications:</h5>
        <p>{data$num_misclass}</p>
        <h5>Observations Misclassified:</h5>
        <pre>{paste(data$which_misclass, collapse = ', ')}</pre>
        <h5>Estimated Mu Values:</h5>
        <pre>{format_mu(data.frame(do.call(rbind, data$mu)), d = length(data$mu)/k)}</pre>
      ")
            box(
              title = paste("Cluster Information for", k, "Clusters"),
              collapsible = TRUE,
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              div(
                style = 'overflow-y: auto; max-height: 400px;',
                HTML(content),
                h5("Plot of Log-Likelihood:"),
                plotOutput(outputId = paste0("loglik_plot_", k))
              )
            )
          })
        )
      })
      
      }else if (grepl("KMM",input$clustering_method,fixed=TRUE)) {
        print("STARTING CLUSTER")
        
        # Perform clustering
        for (cluster in 1:max_clusters) {
          cluster_results[[cluster]] <- kmmeans(as.data.frame(log(clustering_data[[2]]+1/6)), cluster,n.init = 100)
        }
        print("FINISH CLUSTER")
        

        for (cluster in 1:max_clusters) {
          print("hehe")
          cluster_results[[cluster]]$ARI <- mclust::adjustedRandIndex(
            cluster_results[[cluster]]$partition,
            mclust::map(clustering_data[[4]])
          )
        }
        print("pizza")
        optimal_cluster <- KMM_optimalG(cluster_results,clustering_data[[2]])
        optimal_data <- cluster_results[[optimal_cluster]]
        optimal_errors <- mclust::classError(
          optimal_data$partition,
          mclust::map(clustering_data[[4]])
        )
        print("2")
        optimal_mu_est <- optimal_data$centers
        optimal_ARI <- optimal_data$ARI
        
        output$cluster_info_boxes <- renderUI({
          req(cluster_results) 
          
          
          print("brirburur")
          cluster_data <- isolate({
            lapply(1:max_clusters, function(k) {
              ARI <- cluster_results[[k]]$ARI
              misclasss <- mclust::classError(
                cluster_results[[k]]$partition,
                mclust::map(clustering_data[[4]])
              )
              list(
                ARI = ARI,
                num_misclass = length(misclasss$misclassified),
                which_misclass = misclasss$misclassified,
                mu = cluster_results[[k]]$centers
              )
            })
          })
          print("mememem")
          # Generate the UI output
          tagList(
            lapply(1:1, function(fgj){
              data <- cluster_data[[optimal_cluster]]
              content <- glue::glue("
        <h5>Optimal Number of Clusters (by criterion):</h5>
        <p>{optimal_cluster}</p>
        <h5>Adjusted Rand Index (ARI):</h5>
        <p>{data$ARI}</p>
        <h5>Number of Misclassifications:</h5>
        <p>{data$num_misclass}</p>
        <h5>Observations Misclassified:</h5>
        <pre>{paste(data$which_misclass, collapse = ', ')}</pre>
      ")
              box(
                title = paste("Summary of Optimal Clustering Information"),
                collapsible = TRUE,
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                div(
                  style = 'overflow-y: auto; max-height: 400px;',
                  HTML(content)
                )
              )
            }),
            lapply(1:max_clusters, function(k) {
              data <- cluster_data[[k]]
              content <- glue::glue("
        <h5>Adjusted Rand Index (ARI):</h5>
        <p>{data$ARI}</p>
        <h5>Number of Misclassifications:</h5>
        <p>{data$num_misclass}</p>
        <h5>Observations Misclassified:</h5>
        <pre>{paste(data$which_misclass, collapse = ', ')}</pre>
        <h5>Estimated Mu Values:</h5>
        <pre>{format_mu(data$mu, d = length(data$mu)/k)}</pre>
      ")
              box(
                title = paste("Cluster Information for", k, "Clusters"),
                collapsible = TRUE,
                width = NULL,
                solidHeader = TRUE,
                status = "primary",
                div(
                  style = 'overflow-y: auto; max-height: 400px;',
                  HTML(content)
                )
              )
            })
          )
        })
    }
  })
  
  
}
