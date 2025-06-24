

# Step 6: Helper functions for data conversion
# ---------------------------------------------

# Convert R data to Python numpy arrays
r_to_numpy <- function(r_data) {
  if (is.vector(r_data)) {
    return(np$array(r_data))
  } else if (is.matrix(r_data) || is.data.frame(r_data)) {
    return(np$array(as.matrix(r_data)))
  } else {
    stop("Unsupported data type")
  }
}

# Convert Python numpy arrays back to R
numpy_to_r <- function(py_array) {
  return(py_to_r(py_array))
}



evpruning_R <- function(data,tau=0.3) {
  
  # Convert R data to Python format
  py_data <- r_to_numpy(data)
  
  
  # Call the Python function
  result <- eigenvalue_pruning(py_data,tau)
  
  return(numpy_to_r(result))
}



tensor_to_r <- function(tensor) {
  # Convert tensor to numpy, then to R
  numpy_array <- tensor$detach()$numpy()
  r_data <- py_to_r(numpy_array)
  return(r_data)
}


que_R <- function(data,tau=0.3) {
  
  # Convert R data to Python format
  py_data <- r_to_numpy(data)
  
  
  # Call the Python function
  result <- que_mean(py_data,tau)
  
  ##Not sure why there is a tensor here
  return(tensor_to_r(result))
}






get_imp_estimate <- function(X, method="norm.nob"){
  
  imp <- mice(X, method = method, m = 5, printFlag = FALSE)
  all_means <- sapply(1:5, function(i) {
    completed_data <- complete(imp, i)
    colMeans(completed_data, na.rm = TRUE)
  })
  
  # Average across imputations
  normimp <- rowMeans(all_means)
  
}

get_estimate <- function(X, method="MMD"){
  
  if (method=="MMD"){
    res<-mmd_estown(X, model = "multidim.Gaussian.loc", par2 = 1, control=list(method="SGD", burnin=100, nstep=300))$estimator
  }else if (method=="mean"){
    
    res<-colMeans(X, na.rm=T)
    
  }else if (method=="median"){
    
    res<-apply(X,2, median, na.rm=T)
    
  }else if (method=="normimp"){
    
    res<-get_imp_estimate(X, method="norm.nob")
  }else if (method=="meanimp"){
    
    res<-get_imp_estimate(X, method="mean")
    
  } else if (method=="evpruning"){
    ## only look at complete cases
    res<-evpruning_R(X[complete.cases(X), ] )
    
  }else if (method=="que"){
   ## only look at complete cases
    res<-que_R(X[complete.cases(X), ] )
  }
  
  
  return(res)
  
}

generate_markov_missingness <- function(X, p1, p2) {
  # X: n x d data matrix
  # p1: probability of missing if previous dimension was missing (or for first dimension)
  # p2: probability of missing if previous dimension was observed (p2 > p1)
  
  n <- nrow(X)
  d <- ncol(X)
  
  # Create a copy of X to modify
  X_missing <- X
  
  # Create missingness indicator matrix (TRUE = missing, FALSE = observed)
  miss_pattern <- matrix(FALSE, nrow = n, ncol = d)
  
  # First dimension: missing with probability p1 (vectorized)
  miss_pattern[, 1] <- runif(n) < p1
  
  miss_pattern[, 2] <- runif(n) < p1
  
  # Subsequent dimensions follow Markov chain (vectorized)
  if (d > 1) {
    for (j in c(5,9)) {
      # Generate random numbers for all observations at once
      rand_vals <- runif(n)
      
      # Vectorized logic: use p2 if previous was missing, p1 if previous was observed
      probs <- ifelse(miss_pattern[, 1], p2, 0)
      miss_pattern[, j] <- rand_vals < probs
    }
    
    
    for (j in c(6,10)) {
      # Generate random numbers for all observations at once
      rand_vals <- runif(n)
      
      # Vectorized logic: use p2 if previous was missing, p1 if previous was observed
      probs <- ifelse(miss_pattern[, 2], p2, 0)
      miss_pattern[, j] <- rand_vals < probs
    }
    
  }
  
  # Apply missingness to data (vectorized)
  X_missing[miss_pattern] <- NA
  
  return(X_missing)
}


generate_MNAR_missingness <- function(X,alpha=0, beta=rep(1,ncol(X)),nu=5) {
  # X: n x d data matrix
  # p1: probability of missing if previous dimension was missing (or for first dimension)
  # p2: probability of missing if previous dimension was observed (p2 > p1)
  ##Continue here!!
  
  n<-nrow(X)
  
  #prob_obs <- pt(alpha + X%*%beta, df = nu)
  prob_obs <- pnorm(alpha + X%*%beta)
  r <- rbinom(n, 1, prob_obs)
  
  # Create observed Y (set missing to NA)
  X_mis <- X
  X_mis[r == 0, ] <- NA
  
  return(X_mis[r == 1, ])
  
  
}


createdistMNAR<-function(n, d, eps, p1,p2, ...){
  
  args<-list(...)
  
  cont<-rbinom(n=n,size=1, prob=eps)
  
  print(cont)
  
  X<-matrix(rnorm(n=n*d), nrow=n, ncol=d ) 
  
  XMCAR <- generate_markov_missingness(X[cont==0,],p1,p2) 
  XMNAR <- generate_MNAR_missingness(X[cont==1,], ...)
  
  return(rbind(XMCAR, XMNAR))
  
  
}

createdistMCAR<- function(n,d,eps, p1=0.05,p2=0.1, contdist, ...){
  
  ## contamination index
  cont<-rbinom(n=n,size=1, prob=eps)
  
  X<-matrix(rnorm(n=n*d), nrow=n, ncol=d ) 
  
  X[cont==1,]<-drawcontamination(n=sum(cont==1),d=d, contdist=contdist, ...)
  
  ##include MCAR!
  X_mis<-generate_markov_missingness(X, p1, p2)
  
  return(X_mis)
  
}




drawcontamination<-function(contdist,n, d, ...){
  
  args<-list(...)
  
  if (contdist=="rnorm"){
    mu=args$mu
    return(matrix(rnorm(n=n*d, mean=mu), nrow=n, ncol=d ))
    
  } else if (contdist=="dirac"){
    c<-args$c
    return(matrix(c, nrow=n, ncol=d ) )
    
  }
  
  
}


createdistMNARMCAR<-function(n, d, eps, p1,p2, contdist, ...){
  
  args<-list(...)
  
  X<-matrix(rnorm(n=n*d), nrow=n, ncol=d ) 
  
  ## contamination index
  cont1<-rbinom(n=n,size=1, prob=eps)
  cont2<-rbinom(n=n,size=1, prob=eps)
  
  # Full data contamination
  X[cont1==1,]<-drawcontamination(n=sum(cont1==1),d=d, contdist=contdist, ...)
  
  
  #Missingness Contamination
  XMCAR <- generate_markov_missingness(X[cont2==0,],p1,p2) 
  XMNAR <- generate_MNAR_missingness(X[cont2==1,])
  
  return(rbind(XMCAR, XMNAR))
  
  
}




