library(MASS)
library(drf)
source("MMD_minimization.R")



############################################
## First: Trivial unconditional example just to test the function ##
###########################################

set.seed(1)

n<-200
p<-2
X<-matrix(rnorm(p*n, 1), nrow=n)


simulate_P_theta <- function(M, theta){ mvrnorm(M, mu=theta, Sigma=diag(length(theta)))   }

derivative_log_density<-function(theta, Y){sweep(Y,2, theta, "-")  }


res<-MMDestimation(X, M=5, k=NULL, w=NULL, numit=2000, simulate_P_theta, derivative_log_density, theta_init=rep(0,p))

res




############################################
## Second: Estimation of conditional mean ##
###########################################
#### Model: Y = mu(x) + N(0,1) ####

set.seed(2)

n<-2000
beta1<-1
beta2<--1.8
d<-1

# Model Simulation
X<-mvrnorm(n = n, mu=c(0,0), Sigma=matrix(c(1,0,0,1), nrow=2,ncol=2))
#X<-mvrnorm(n = n, mu=c(0,0), Sigma=diag(2))
u<-rnorm(n=n, sd = 1)#sqrt(exp(X[,1]))
Y<- matrix(beta1*X[,1] + beta2*X[,2] + u, ncol=1)

# Choose an x that is not too far out
x<-matrix(c(1,1),ncol=2)





## Fit the new DRF framework
drf_fit <- drf(X=X, Y=Y, num.trees=2000, min.node.size = 15,  ci.group.size = 1)

## predict weights
DRF = predict(drf_fit, newdata=x)
weights <- DRF$weights[1,]

#sum(weights*Y)


simulate_P_theta <- function(M, theta){ mvrnorm(M, mu=theta, Sigma=diag(length(theta)))   }

derivative_log_density<-function(theta, Y){sweep(Y,2, theta, "-")  }


res<-MMDestimation(Y, M=5, k=NULL, w=weights, numit=1000, simulate_P_theta, derivative_log_density, theta_init=rep(0,d))


paste("Estimate: ", res)
paste("Truth: ", beta1*x[,1] + beta2*x[,2])

# ############################################
# ## Adding Confidence intervals (takes a bit longer to run because of bad implementation) ##
# ###########################################
# ## Fit the new DRF framework
# drf_fit <- drfCI(X=X, Y=Y, num.trees=500, B=50)
# 
# ## predict weights
# DRF = predictdrf(drf_fit, x=x)
# 
# #sum(weights*Y)
# 
# 
# simulate_P_theta <- function(M, theta){ mvrnorm(M, mu=theta, Sigma=diag(length(theta)))   }
# 
# derivative_log_density<-function(theta, Y){sweep(Y,2, theta, "-")  }
# 
# 
# res<-MMDestimation(Y, M=5, k=NULL, w=DRF$weights[1,], numit=200, simulate_P_theta, derivative_log_density, theta_init=rep(0,d))
# res
# 
# CIdat<-sapply(1: length(DRF$weightsb), function(b) 
#   
#   MMDestimation(Y, M=5, k=NULL, w=DRF$weightsb[[b]][1,], numit=200, simulate_P_theta, derivative_log_density, theta_init=rep(0,d))
#   
#   )
# 
# CIupper <- res+qnorm(1-0.05/2)*sqrt(var(CIdat))
# CIlower <- res-qnorm(1-0.05/2)*sqrt(var(CIdat))
# 
# 
# c(CIlower, res, CIupper)


############################################
## Third: Estimation of conditional beta ##
###########################################

#### Model: Y_2 = beta(x)* Y_1 + N(0,1) ####
#### (Y_1,Y_2)|X=x ~ N(0, [1, beta(x)\\ beta(x), beta(x)^2 + 1  ])


###Double check this!!

set.seed(2)

n<-2000
beta1<-1
beta2<--1.8
d<-1

# Model Simulation
X<-mvrnorm(n = n, mu=c(0,0), Sigma=matrix(c(1,0,0,1), nrow=2,ncol=2))
Y<-matrix(NaN, nrow=n,ncol=2)
betaX<- beta1*X[,1] + beta2*X[,2]

Y<-t ( sapply(1:n, function(j){
mvrnorm(1, mu=c(0,0), Sigma= matrix( c(1, betaX[j], betaX[j], betaX[j]^2 + 1 ), nrow=2,ncol=2,byrow = T  ) )
}) )
# Choose an x that is not too far out
x<-matrix(c(1,1),ncol=2)




## Fit the new DRF framework
drf_fit <- drf(X=X, Y=Y, num.trees=2000, min.node.size = 5,  ci.group.size = 1)

## predict weights
DRF = predict(drf_fit, newdata=x)
weights <- DRF$weights[1,]

#sum(weights*Y)


simulate_P_theta <- function(M, theta){ mvrnorm(M, mu=c(0,0), Sigma= matrix( c(1, theta, theta, theta^2 + 1 ), nrow=2,ncol=2,byrow = T  ) )   }


derivative_log_density<-function(theta, Y){
  
  # Inversion formula: #1/(theta^2 + 1 - theta^2) *c(theta^2 + 1, -theta, -theta, 1 )
  Siginv<-matrix( c(theta^2 + 1, -theta, -theta, 1 ), nrow=2,ncol=2,byrow = T  )
  Derivmatrix<-matrix( c(0, 1, 1, 2*theta), nrow=2,ncol=2,byrow = T  )
  # d L/ dSigma * dSigma/dtheta
  return( 
    matrix(sapply(1:nrow(Y), function(j) -0.5*sum(diag( t(Siginv  - Siginv%*%Y[j,]%*%t(Y[j,])%*%Siginv )%*%Derivmatrix)) ), ncol=1)
  )
}


res<-MMDestimation(Y, M=10, k=NULL, w=weights, numit=1000, simulate_P_theta, derivative_log_density, theta_init=0)


paste("Estimate: ", res)
paste("Truth: ", beta1*x[,1] + beta2*x[,2])


############################################
## Fourth: Estimation of a unique beta but conditional variance! ##
###########################################
###Continue here!!
#### Model: Y_2 = beta Y_1 + sigma(x)*N(0,1) ####
#### (Y_1,Y_2)|X=x ~ N(0, [1, beta\\ beta, beta^2 + sigma(x)^2  ])





set.seed(2)

n<-1000
beta1<-1
beta2<--1.8
d<-1

# Model Simulation
X<-mvrnorm(n = n, mu=c(0,0), Sigma=matrix(c(1,0,0,1), nrow=2,ncol=2))
#X<-mvrnorm(n = n, mu=c(0,0), Sigma=diag(2))
beta<- c(1,2)

Y<-t ( sapply(1:n, function(j){
  mvrnorm(1, mu=c(0,0), Sigma= matrix( c(1, beta, beta, beta^2 + exp(-0.2*X[,1]) ), nrow=2,ncol=2,byrow = T  ) )
}) )

# Choose an x that is not too far out
x<-matrix(c(1,1),ncol=2)




## Fit the new DRF framework
drf_fit <- drf(X=X, Y=Y, num.trees=2000, min.node.size = 15,  ci.group.size = 1)

## predict weights
DRF = predict(drf_fit, newdata=X)
weights <- DRF$weights

#sum(weights*Y)

####Continue here!!!
simulate_P_theta <- function(M, theta, nuisance){ 
  
  # Step 1: Bootstrap from nuisance parameters (trying to make it more general)
  sigx2=nuisance[sample(1:n, M, replace=T),]
  
  # Step 2: Simulate from Y|X_i
  return(mvrnorm(M, mu=c(0,0), Sigma= matrix( c(1, theta, theta, theta^2 + sigx2 ), nrow=2,ncol=2,byrow = T  ) ) )   }



derivative_log_density<-function(theta, Y){
  
  # Inversion formula: #1/(theta^2 + 1 - theta^2) *c(theta^2 + 1, -theta, -theta, 1 )
  Siginv<-matrix( c(theta^2 + 1, -theta, -theta, 1 ), nrow=2,ncol=2,byrow = T  )
  Derivmatrix<-matrix( c(0, 1, 1, 2*theta), nrow=2,ncol=2,byrow = T  )
  # d L/ dSigma * dSigma/dtheta
  return( 
    matrix(sapply(1:nrow(Y), function(j) -0.5*sum(diag( (Siginv  - Siginv%*%Y[j,]%*%t(Y[j,])%*%Siginv )%*%Derivmatrix)) ), ncol=1)
  )
}


res<-MMDestimation(Y, M=5, k=NULL, w=weights, numit=1000, simulate_P_theta, derivative_log_density, theta_init=0)


paste("Estimate: ", res)
paste("Truth: ", sin(beta1*x[,1] + beta2*x[,2]))










