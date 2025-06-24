# Required library
library(kernlab)



# Algorithm 1
MMDestimation <- function(X, M=10, k=NULL, w=NULL, numit, simulate_P_theta, derivative_log_density, theta_init, nuisance=NULL) {
  
  
  # nuisance: matrix of nrow(Y) \times # of nuisance parameter that can change with X. Default=NULL
  # simulate_P_theta(M, theta): Inputs M and theta, simulates M i.i.d. realization from P_{theta}
  
  # derivative_log_density(theta, Y): Inputs theta and Y and evaluates the 
  # derivative of the log density. Returns a matrix with each row j representing the derivative of the log density for Y_j.
  # Example for the Gaussian mean: sapply(1:M, function(j) solve(Sigma)*(Y[j,]-mu))
  
  n<-nrow(X)
  d<-ncol(X)
  
  
  
  if (is.null(w)){
    #unconditional estimation
    w <- rep(1/n, n)
  }
  
  
  theta<-matrix(NA, nrow=numit+1, length(theta_init))
  # Initialize theta
  theta[1,] <- theta_init
  
  # Kernel function
  if (is.null(k)){
    k <-rbfdot(sigma=1)
  }
  
  ### Check how to best choose eta!! (page 17 of MMD estimation papers)
  eta<-p^2*M/sqrt(numit) ##p^2 seems to improve!
  
  # SGD loop
  #while( dist(rbind(thetaold, theta)) > eps ) {
   for (t in 1:numit){ 
     
     if (t%%500==0){
     cat(t)
     }
     ###Choose eta!!####
     
     # Generate samples
     if (is.null(nuisance)){
     Y <-  simulate_P_theta(M,theta[t,])
     }else{
     Y <-  simulate_P_theta(M,theta[t,], nuisance)
     }
     
     ##Implement things over patterns!!!

     KYY<-kernelMatrix(k, x=Y, y = Y)
     KXY<-kernelMatrix(k, x=X, y = Y)
     
     KYYnodiag<-KYY
     diag(KYYnodiag)<-0

     meanY<-1/(M-1)*rowSums(KYYnodiag)#apply(KYYnodiag,1, sum)
     meanXY<-colSums(w*KXY) #for w=(1/n,...,1/n) this should just be colMeans
    
    # Compute gradient
     if (is.null(nuisance)){
     gradientY <- derivative_log_density(theta[t,], Y)
     }else{
     gradientY <- derivative_log_density(theta[t,], Y, nuisance) 
     }
    # Combine everything
    Crit<-2/M*colSums((meanY - meanXY)*gradientY)
    ###
  
    #thetaold<-theta
    
    # Update theta
    theta[t+1,] <- theta[t,] - eta * Crit
  }
  

  # Return the final estimate, maybe weight more recent iterations more
  #return(colMeans(theta[2:(numit+1),, drop=F]))
  return(colMeans(theta[2:(numit+1),, drop=F]))
  
  
}

drfCI <- function(X, Y, B, sampling = "binomial",...) {
  
  n <- dim(X)[1]
  
  # compute point estimator and DRF per halfsample
  # weightsb: B times n matrix of weights
  DRFlist <- lapply(seq_len(B), function(b) {
    
    # half-sample index
    indexb <- if (sampling == "binomial") {
      seq_len(n)[as.logical(rbinom(n, size = 1, prob = 0.5))]
    } else {
      sample(seq_len(n), floor(n / 2), replace = FALSE)
    }
    
    ## Using normal Bootstrap on the data and refitting DRF
    DRFb <- 
      drf(X = X[indexb, , drop = F], Y = Y[indexb, , drop = F],
          ci.group.size = 1, ...)
    
    
    return(list(DRF = DRFb, indices = indexb))
  })
  
  return(list(DRFlist = DRFlist, X = X, Y = Y) )
}


predictdrf<- function(DRF, x, functional = NULL, ...) {
  
  require(Matrix)
  
  ntest <- nrow(x)
  n <- nrow(DRF$Y)
  
  
  weightsb <- lapply(DRF$DRFlist, function(l) {
    
    weightsbfinal <- Matrix(0, nrow = ntest, ncol = n , sparse = TRUE)
    
    weightsbfinal[, l$indices] <- predict(l$DRF, x)$weights 
    
    return(weightsbfinal)
  })
  
  
  weightsall <- Reduce("+", weightsb) / length(weightsb)
  
  if (!is.null(functional)) {
    ## Achtung: This so far works only for one x!  
    stopifnot("Not yet implemented for several x" = ntest == 1)
    
    thetahatb <- 
      lapply(weightsb, function(w)  
        functional(weights = w , X = DRF$X, Y = DRF$Y, x = x))
    thetahatbforvar <- do.call(rbind, thetahatb)
    thetahat <- functional(weights = weightsall , X = DRF$X, Y = DRF$Y, x = x)
    thetahat <- matrix(thetahat, nrow = dim(x)[1])
    var_est <- if (dim(thetahat)[2] > 1){  
      a <- sweep(thetahatbforvar, 2, thetahat, FUN = "-")
      crossprod(a, a) / B
    } else {
      mean((c(thetahatbforvar) - c(thetahat)) ^ 2) 
    }
    
    return(list(weights = weightsall, thetahat = thetahat, weightsb = weightsb, 
                var_est = var_est ))
    
  } else {
    return(list(weights = weightsall, weightsb = weightsb ))
  }
}


















