##Functions 

generate_data <- function(par,N){
  names(par)<-c("pi","mu","Sigma")
  
  d <- length(par[[2]][[1]])
  
  Y<-matrix(0,ncol=d,nrow=N)
  X<-matrix(0,ncol=d,nrow=N)
  
  
  z<-t(rmultinom(n=N,size=1,prob=par[[1]]))
  
  for (i in 1:N){
    grp<-which(z[i,]==1)
    X[i,]<-rmvnorm(1,mean=par[[2]][[grp]],sigma=par[[3]][[grp]])  
    for (j in 1:d){
      Y[i,j]<-rpois(1,exp(X[i,j]))
    }
  }
  
  dat<-list()
  dat[[1]]<-par
  dat[[2]]<-Y
  dat[[3]]<-X
  dat[[4]]<-z
  
  return(dat)
}

mixPLN <-function(dat,G){
  ptm<-proc.time()
  ###### Parameter Updates ####
  step<- 0.001
  ###Reading data from simulation study. Y is the observed counts. true contains the true group label and true_par is the true values of the parameters that we need to assess performance later and 
  Y<-dat[[2]]
  true<-dat[[4]]
  true_par<-dat[[1]]
  
  #true_par
  
  
  N<-nrow(Y) #sample size
  d<-ncol(Y) #dimension of the data
  
  
  O_mat <- matrix(NA,nrow = N,ncol = d)
  for(col in 1:d){
    for (row in 1:N){
      if(!is.na(Y[row,col])){
        O_mat[row,col] <- 1
      }
    }
  }
  
  #This should be redundant due to the way pattern was created but leaving here in case
  all_miss <- apply(O_mat, 1, function(x) all(is.na(x)))
  O_mat<- O_mat[ !all_miss, ]
  

  #updated count matrix with entries determined to be missing set to 0.
  #aggr(O_mat, plot = TRUE, numbers = TRUE, prop = FALSE)
  Y[is.na(O_mat)]<-NA
  
  ###Removes any observation with all values missing
  all_miss <- apply(Y, 1, function(x) all(is.na(x)))
  Y<- Y[ !all_miss, ]
  
  N<-nrow(Y) #sample size
  d<-ncol(Y) #dimension of the data
  
  O_list<-list()
  for (i in 1:N){
    O_list[[i]]<-diag(d)[which(O_mat[i,]==1),,drop=FALSE]
  }
  
  Y[is.na(Y)]<--999
  
  all_miss <- apply(Y, 1, function(x) all(is.na(x)))
  Y<- Y[ !all_miss, ]
  
  lib_mat<-rep(1,d) #for simulation, we leave the lib_mat as a vector of 1. but for real data, we uncomment the one below.
  
  #lib_mat<-edgeR::calcNormFactors(Y)

  #### Initialization ###
  mu<-list()
  psi<-list()
  lambda<-list()
  sigma<-list()
  isigma<-list()
  sigma_new<-list()
  m<-list()
  S<-list()
  P<-list()
  Q<-list()
  
  ###Other intermediate items initialized
  start<-list()
  Sk<-array(0, c(d,d,G) )
  GX<-list()
  dGX<-list()
  iOsigO<-list()
  
  z_S<-list()
  z_SO<-list()
  z_DO<-list()
  
  ini_Y<-Y
  ini_Y[ini_Y==-999]<-NA
  rownames(ini_Y)<-1:nrow(ini_Y)
  kept<-as.numeric(rownames(na.omit(ini_Y)))
  k_means<-kmmeans(as.data.frame(log(ini_Y+1)),K=G,n.init = 100)[["partition"]]  ##Using k missing means to start to algorithm
  z<-mclust::unmap(k_means) ##Initial value for Z
  pi_g<-colSums(z)/N
  
  ###Initial value for Mu and Sigma
  ini_Y2<-na.omit(ini_Y)
  ini_Y[is.na(ini_Y)]<-0
  
  ###Initial value for Mu and Sigma
  for (g in 1:G){
    obs<-which(z[kept,g]==1)
    mu[[g]]<-colMeans(log(ini_Y2[obs,]+1/6)) 
    sigma[[g]]<-var(log(ini_Y2[obs,]+1/6))
    isigma[[g]]<-solve(sigma[[g]])
  }
  

  ###Initial value for m and S
  for (g in 1:G){
    S[[g]]<-list()
    start[[g]]<-list()
    m[[g]]<-list()
    iOsigO[[g]]<-list()
    for (i in 1:N){
      do<-nrow(O_list[[i]])
      start[[g]][[i]]<-log(O_list[[i]]%*%ini_Y[i,]+1/6) ###Starting value for M
      m[[g]][[i]]<-start[[g]][[i]]
      S[[g]][[i]]<-diag(do)*0.000000001
      iOsigO[[g]][[i]]<-solve(as.matrix(O_list[[i]])%*%as.matrix(sigma[[g]])%*%t(as.matrix(O_list[[i]])))
    }
  }
  
  checks<-0
  it<-1
  aloglik<-NULL
  loglik<-NULL ##Log likelihood is stored in this vector to check for convergence.
  aloglik[c(1,2,3)]<-0
  
  for(g in 1:G){
    start[[g]] <- lapply(start[[g]],as.vector)
    m[[g]] <- lapply(m[[g]],as.vector)
    
  }
  
  
  EM_list<-EM_loop(start, lib_mat, m, mu, sigma, isigma, iOsigO, S, O_list, 
                   as.matrix(Y), as.matrix(z), as.matrix(O_mat), pi_g, d, N, G)
  
  loglik <- as.vector(EM_list[["loglik"]])
  z <- as.matrix(EM_list[["z"]])
  pi_g <- unlist(EM_list[["pi_g"]])
  mu <- EM_list[["mu"]]
  sigma <- EM_list[["sigma"]]
  exit_code <- EM_list[["exit_code"]]
  
  
  
  k<- (G-1) +G*d +G*d*(d+1)/2
  BIC <- 2*loglik[length(loglik)] - k*log(N)   
  #plot(loglik,type="l")
  ptm<-proc.time() - ptm
  Y[which(Y==-999)]<-NA
  return(list(Y=Y,pi_g=pi_g,mu=mu,sigma=sigma,z=z,loglik=loglik,BIC=BIC,kmeans=k_means,true=true,time=ptm))
}

generate_missing <- function(dat,percent_miss = 0.10,mech="MAR"){
  local_data <- dat
  
  if(grepl("MAR",mech,fixed=TRUE)){
    mech <- "MAR"
  }else if(grepl("MNAR",mech,fixed=TRUE)){
    mech<-"MNAR"
  }else if(grepl("MCAR",mech,fixed=TRUE)){
    mech<-"MCAR"
  }else{
    print("Unrecognized mech inpuit, deafult is MAR")
    mech<-"MAR"
  }
  print("PASSED MECH")
  d<- ncol(local_data[[2]])
  local_data[[5]] <- local_data[[2]] #index 5 is for complete data, index 2 is to be made partially missing
  #percent_miss functions slightly differently now. It is the percent of observations which are incomplete
  #as a rough rule, the number of missing cells will be about 0.4 times percent_miss
  pattern <- matrix(1,ncol=d)
  for(i in seq(d,2,-1)){
    pattern <- rbind(pattern,gtools::permutations(d,d,c(rep(1,times=i),rep(0,times=d-i)),set=FALSE))
  }
  pattern <- distinct(as.data.frame(pattern))
  pattern <- pattern[-1,]
  print("MADE PATTERN")
  #Sometimes ampute can generate columns of all missing values. This will cause method to break.
  #We implement check, and if this occurs we just regenerate missing values
  #Rare enough where this shouldn't impact performance drastically
  attempts <- 1
  while(attempts < 10){
    candidate_missing_dat <- as.matrix(ampute(local_data[[2]],prop=percent_miss,patterns = pattern,mech = mech)[["amp"]])
    for(col in 1:ncol(candidate_missing_dat)){
      if(all(is.na(candidate_missing_dat[,col]))){
        #we have found col with all missing values, restart loop without updating run
        attempts <- attempts + 1
        next
      }
    }
    #If we loop through all cols and there are no all missing cols we are good
    local_data[[2]] <- candidate_missing_dat
    print("DATA AMPUTED")
    break
  }
  if(attempts >= 10){
    stop("Missing data process is consistently creating columns of all missing values!")
  }
  
  return(local_data)
}

format_params <- function(parameters){
  alphas <- c()
  sd_lowers <- c()
  sd_uppers <- c()
  
  mus <- list()
  pi_g <- c()
  
  for(cluster in 1:length(parameters)){
    alphas[cluster] <- parameters[[cluster]]$alpha
    sd_lowers[cluster] <- parameters[[cluster]]$sd_lower
    sd_uppers[cluster] <- parameters[[cluster]]$sd_upper
    pi_g[cluster] <- parameters[[cluster]]$percent
    
    mus[[cluster]] <- parameters[[cluster]]$mu
  }
  
  return_par <- list()
  
  return_par[[1]] <- pi_g
  return_par[[2]] <- mus
  return_par[[3]] <- generate_sigmas(parameters)
  
  return(return_par)
  
}

generate_sigmas <- function(params){
  return_list <- list()
  d<- length(params[[1]]$mu)
  for (i in seq_along(params)) {
    stop<-0
    while (stop<1){
      cor <- rcorrmatrix(d,params[[i]]$alpha)
      sd <- runif(d,params[[i]]$sd_lower,params[[i]]$sd_upper)
      sigma <- diag(sd)  %*% cor  %*% diag(sd)
      if (min(eigen(sigma)$values)>0.3) stop<-1
    }
    return_list[[i]] <- sigma
  }
  return(return_list)
}


plot_loglikelihood <- function(loglik_values, num_clusters) {
  # Convert log-likelihood values to a dataframe for ggplot
  loglik_df <- data.frame(
    Iteration = seq_along(loglik_values),  
    LogLikelihood = loglik_values       
  )
  
  # Create the ggplot
  my_plot <- ggplot(loglik_df, aes(x = Iteration, y = LogLikelihood)) +
    geom_line(color = "blue", size = 1) +  
    geom_point(color = "red", size = 2) +  
    labs(
      title = paste("Log-Likelihood for", num_clusters, "Clusters"),
      x = "Iteration",
      y = "Log-Likelihood"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5)  
    )
  
  return(my_plot)
}


KMM_optimalG <- function(cluster_results,log_count_mat){
  n<- dim(log_count_mat)[1]
  p_bar <- sum(!is.na(log_count_mat))/n
  W_vector <- vector()
  for(clusters in 1:length(cluster_results)){
    W_vector[clusters] <- cluster_results[[clusters]]$criterion
  }
  D_vector <- vector()
  D_vector[1] <- 0
  for(clusters in 1:length(cluster_results)){
    D_vector[clusters+1] <- W_vector[clusters]/(n*p_bar)
  }
  J_vector <- vector()
  for(clusters in 1:length(cluster_results)){
    J_vector[clusters] <- (D_vector[clusters+1])^(-p_bar/2) - (D_vector[clusters])^(-p_bar/2)
  }
  return(which.max(J_vector))
}

#ImNotaGit on stackoverflow
scat.my <- function(data, mapping, colors, ...) {
  # Parse x and y variable names
  x <- as.character(unclass(mapping$x))[2]
  y <- as.character(unclass(mapping$y))[2]
  
  # Extract and filter data
  plottingdat1 <- data.table(x = data[[x]], y = data[[y]], Cluster = data$Cluster)[x != -666 & y != -666]
  
  ggplot(plottingdat1, aes(x = x, y = y, color = as.factor(Cluster))) +
    geom_point() +
    scale_color_manual(values = colors) + 
    theme_minimal()
}

diag.my <- function(data, mapping, colors, ...) {
  x <- as.character(unclass(mapping$x))[2]
  
  plottingdat2 <- data.table(x = data[[x]], Cluster = data$Cluster)[x != -666]
  
  # Check if data is empty after filtering
  if (nrow(plottingdat2) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No data to plot", size = 6) + theme_void())
  }
  
  # Create density plot with cluster-based coloring
  ggplot(plottingdat2, aes(x = x, color = as.factor(Cluster))) +
    geom_density(alpha = 0.7) +
    scale_color_manual(values = colors) +
    theme_minimal()
}

format_mu <- function(mu_matrix, d) {
  # Flatten the matrix into a vector
  mu_vector <- as.vector(mu_matrix)  
  
  
  num_values <- length(mu_vector)
  num_rows <- ceiling(num_values / d)  
  padded_vector <- c(mu_vector, rep("", d * num_rows - num_values))  
  matrix_form <- matrix(padded_vector, nrow = num_rows, byrow = TRUE)
  
  
  formatted <- apply(matrix_form, 1, function(row) paste(row, collapse = "\t"))
  paste(formatted, collapse = "\n")
}
