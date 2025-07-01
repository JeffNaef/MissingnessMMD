

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

# generate_markov_missingness <- function(X, p1=0.2, p2=1, p3=0.1) {
#   # X: n x d data matrix
#   # p1: probability of missing if previous dimension was missing (or for first dimension)
#   # p2: probability of missing if previous dimension was observed (p2 > p1)
#   
#   n <- nrow(X)
#   d <- ncol(X)
#   
#   # Create a copy of X to modify
#   X_missing <- X
#   
#   # Create missingness indicator matrix (TRUE = missing, FALSE = observed)
#   miss_pattern <- matrix(FALSE, nrow = n, ncol = d)
#   
#   # First dimension: missing with probability p1 (vectorized)
#   miss_pattern[, 1] <- runif(n) < p1
#   miss_pattern[, 2] <- runif(n) < p1
#   miss_pattern[, 3] <- runif(n) < p3
#   
# 
#   
#   # Subsequent dimensions follow Markov chain (vectorized)
#   if (d > 1) {
#     for (j in c(5,9)) {
#       # Generate random numbers for all observations at once
#       rand_vals <- runif(n)
#       
#       # Vectorized logic: use p2 if previous was missing, p1 if previous was observed
#       probs <- ifelse(miss_pattern[, 1], p2, 0)
#       miss_pattern[, j] <- rand_vals < probs
#     }
#     
#     
#     for (j in c(6,10)) {
#       # Generate random numbers for all observations at once
#       rand_vals <- runif(n)
#       
#       # Vectorized logic: use p2 if previous was missing, p1 if previous was observed
#       probs <- ifelse(miss_pattern[, 2], p2, 0)
#       miss_pattern[, j] <- rand_vals < probs
#     }
#     
#   }
#   
#   
#   
#   # Apply missingness to data (vectorized)
#   X_missing[miss_pattern] <- NA
#   
#   return(X_missing)
# }
# 
# 
# generate_MNAR_missingness <- function(X,alpha=0, beta=rep(1,ncol(X)),nu=5) {
#   # X: n x d data matrix
#   # p1: probability of missing if previous dimension was missing (or for first dimension)
#   # p2: probability of missing if previous dimension was observed (p2 > p1)
#   ##Continue here!!
#   
#   n<-nrow(X)
#   
#   #prob_obs <- pt(alpha + X%*%beta, df = nu)
#   prob_obs <- pnorm(alpha + X%*%beta)
#   r <- rbinom(n, 1, prob_obs)
#   
#   # Create observed Y (set missing to NA)
#   X_mis <- X
#   X_mis[r == 0, ] <- NA
#   
#   return(X_mis[r == 1, ])
#   
#   
# }



generate_markov_missingness <- function(X) {
  # X: n x d data matrix
  # p1: probability of missing if previous dimension was missing (or for first dimension)
  # p2: probability of missing if previous dimension was observed (p2 > p1)
  
  n <- nrow(X)
  d <- ncol(X)
  
  if (!d==10){
    stop("Code only made for d=10 dimensions")
  }
  
  # Create a copy of X to modify
  X_missing <- X
  
  # Create missingness indicator matrix (TRUE = missing, FALSE = observed)
  miss_pattern <- matrix(0, nrow = n, ncol = d)
  
  # First dimension: missing with probability p1 (vectorized)
  
  m<-sample(c(1,2,3,4), size=n, replace = T, prob = c(1/4,1/4,1/4,1/4))
  
  miss_pattern[m==1,]<-matrix(c(1,1,1, rep(0,d-3)), nrow=sum(m==1), ncol=d, byrow = T)
  miss_pattern[m==2,]<-matrix(c(0,0,0, 1,1,1, rep(0,d-6)) , nrow=sum(m==2), ncol=d, byrow = T)
  miss_pattern[m==3,]<-matrix(c(rep(0,d-4),1,1,1,0), nrow=sum(m==3), ncol=d, byrow = T)
  miss_pattern[m==4,]<-matrix(c(rep(0,d)), nrow=sum(m==4), ncol=d, byrow = T)
  
  
  # Apply missingness to data (vectorized)
  X_missing[miss_pattern==1] <- NA
  
  return(X_missing)
}


generate_MNAR_missingness <- function(X,alpha=0, beta=rep(1,ncol(X))/sqrt(ncol(X)), eps=0, contdist=NULL, ...) {
  # X: n x d data matrix
  # p1: probability of missing if previous dimension was missing (or for first dimension)
  # p2: probability of missing if previous dimension was observed (p2 > p1)
  ##Continue here: Need to replace the pnorm with the actual mixture distribution!
  ## For instance in the case of of N(\mu, I) contamination:
  ## (1-eps)*N(0,I) + eps*N(\mu, I)
  
  n<-nrow(X)
  d<-ncol(X)
  args<-list(...)
  
  X_missing<-X
  
  
  if (contdist=="rnorm"){
  
    cdffunc<-function(x, args){
      mu=args$mu
      # X ~ N(mu, I), beta=rep(1,ncol(X))/sqrt(ncol(X)), beta'*X ~ N(mu%*%beta, beta%*%beta)
      return( (1-eps)*pnorm(alpha + x%*%beta, sd=sqrt(beta%*%beta)) + eps*pnorm(alpha + x%*%beta, mean=alpha + rep(mu,d)%*%beta, sd=sqrt(beta%*%beta))  )
    }
  
    
  } else if (contdist=="dirac"){
    
    
    cdffunc<-function(x, args){
      c<-args$c
      # X ~ dirac_c, beta%*%x ~ dirac_beta%*%c 
      return( (1-eps)*pnorm(alpha + x%*%beta, sd=sqrt(beta%*%beta)) + eps*(alpha + x%*%beta >= rep(c,d)%*%beta)  )
    }
  
    
  }else if (is.null(contdist)){
    
    cdffunc<-function(x, args=NULL){
      return(pnorm(alpha + x%*%beta)/2)
    }
  }
  
  
  m <-apply(X,1, function(x) {
      #probs <- c(pnorm(alpha + x%*%beta)/2, 1/2-pnorm(alpha + x%*%beta)/2, 1/2- pnorm(alpha + x%*%beta)/2, pnorm(alpha + x%*%beta)/2)
    probs <- c(cdffunc(x, args)/2, 1/2-cdffunc(x, args)/2, 1/2- cdffunc(x, args)/2, cdffunc(x, args)/2)  
    return(sample(1:4, 1, prob = probs))
    })
  
  # Create missingness indicator matrix (TRUE = missing, FALSE = observed)
  miss_pattern <- matrix(TRUE, nrow = n, ncol = d)
  
  
  miss_pattern[m==1,]<-matrix(c(1,1,1, rep(0,d-3)), nrow=sum(m==1), ncol=d, byrow = T)
  miss_pattern[m==2,]<-matrix(c(0,0,0, 1,1,1, rep(0,d-6)) , nrow=sum(m==2), ncol=d, byrow = T)
  miss_pattern[m==3,]<-matrix(c(rep(0,d-4),1,1,1,0), nrow=sum(m==3), ncol=d, byrow = T)
  miss_pattern[m==4,]<-matrix(c(rep(0,d)), nrow=sum(m==4), ncol=d, byrow = T)
  
  X_missing[miss_pattern==1] <- NA
  
  return(X_missing)
  
  
}


createdistMNAR<-function(n, d, eps, ...){
  
  args<-list(...)
  
  ##MNAR contamination
  cont<-rbinom(n=n,size=1, prob=eps)
  
  X<-matrix(rnorm(n=n*d), nrow=n, ncol=d ) 
  
  X[cont==0,] <- generate_markov_missingness(X[cont==0,]) #MCAR case
  X[cont==1,] <- generate_MNAR_missingness(X[cont==1,],eps=0, contdist=NULL, ...) #MNAR case
  
  return(X)
  
}

createdistMCAR<- function(n,d,eps, contdist, ...){
  
  ## data contamination index
  cont<-rbinom(n=n,size=1, prob=eps)
  
  X<-matrix(rnorm(n=n*d), nrow=n, ncol=d ) 
  
  X[cont==1,]<-drawcontamination(n=sum(cont==1),d=d, contdist=contdist, ...)
  
  ##include MCAR!
  X_mis<-generate_markov_missingness(X)
  
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


createdistMNARMCAR<-function(n, d, eps, contdist, ...){
  
  args<-list(...)
  
  X<-matrix(rnorm(n=n*d), nrow=n, ncol=d ) 
  
  #Independent contamination!
  
  ## data contamination index
  cont1<-rbinom(n=n,size=1, prob=eps)
  ## MCAR contamination
  cont2<-rbinom(n=n,size=1, prob=eps)
  
  # Data contamination
  X[cont1==1,]<-drawcontamination(n=sum(cont1==1),d=d, contdist=contdist, ...)
  
  # Missingness contamination
  X[cont2==0,] <- generate_markov_missingness(X[cont2==0,]) #MCAR case
  X[cont2==1,] <- generate_MNAR_missingness(X[cont2==1,],eps=eps, contdist=contdist, ...)#MNAR case
  
  return(X)
  
  
  
}


# 
# subtractive_noise <- function(data, eps, true_mean=0) {
#   # Subtractive Noise
#   # Note: This follows a different interface because it is not additive
#   
#   n <- nrow(data)
#   d <- ncol(data)
#   
#   temp <- rnorm(d)
#   v <- temp / norm(temp, type = "2")
#   
#   # Project each datapoint onto v through broadcasting
#   projected_data <- as.vector((data - matrix(true_mean, nrow = n, ncol = d, byrow = TRUE)) %*% v)
#   
#   sorted_indices <- order(projected_data)
#   projected_data <- projected_data[sorted_indices]
#   data <- data[sorted_indices, ]  # Store data based on magnitude of projection
#   
#   return(data[1:ceiling((1 - eps) * n), ])
# }



