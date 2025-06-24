

introduce_missing_values <- function(X, p) {
  # Check input validity
  if (!is.matrix(X)) {
    stop("Input must be a matrix")
  }
  
  if (p < 0 || p > 1) {
    stop("Probability p must be between 0 and 1")
  }
  
  # Create a matrix of the same dimensions as X with logical values
  # TRUE indicates the value will be replaced with NA
  missing_mask <- matrix(
    runif(nrow(X) * ncol(X)) < p, 
    nrow = nrow(X), 
    ncol = ncol(X)
  )
  
  # Create a copy of X to avoid modifying the original
  X_with_missing <- X
  
  # Replace values with NA where missing_mask is TRUE
  X_with_missing[missing_mask] <- NA
  
  return(X_with_missing)
}



library(MASS)
library(drf)
source("MMD_minimization.R")



############################################
## First: Trivial unconditional example just to test the function ##
###########################################

set.seed(1)

n<-200
p<-4
X<-matrix(rnorm(p*n, 1), nrow=n)


simulate_P_theta <- function(M, theta){ mvrnorm(M, mu=theta, Sigma=diag(length(theta)))   }

derivative_log_density<-function(theta, Y){sweep(Y,2, theta, "-")  }

###Something is not right for p > 10!!

##Without Missing Values
res<-MMDestimation(X, M=100, k=NULL, w=NULL, numit=1000, simulate_P_theta, derivative_log_density, theta_init=rep(0,p))

res



## Now: Introduce MCAR missing values

X.NA<-introduce_missing_values(X,p=0.3)

##Impute badly
X.imp<-X.NA
X.imp[is.na(X.imp)]<-10

res<-MMDestimation(X.imp, M=100, k=NULL, w=NULL, numit=1000, simulate_P_theta, derivative_log_density, theta_init=rep(0,p))

res




MMDestimation(X.NA, M=100, k=NULL, w=NULL, numit=1000, simulate_P_theta, derivative_log_density, theta_init=rep(0,p))









##Trying around with univariate Gaussian Example
n<-10000
X<-rnorm(n=n)

# P(X_1 missing| X_1)= 1\{X_1 > 2 \}*0.2 + 1\{X_1 <= 2 \}*0.8 
# P(X_1 observed | X_1) = 1\{X_1 > 2 \}*0.8 + 1\{X_1 <= 2 \}*0.2
Missinggivenx1<-function(X){ 
  
  M<-rep(0, length(X))
  M[X>2] <- rbinom(n=length(M[X>2]), size=1, prob=0.2)
  M[X<=2] <- rbinom(n=length(M[X<=2]), size=1, prob=0.8)
  
return(M)
  
  }


# P(X_1 observed | X_1) = exp(X_1)/(1+exp(X_1)) => the smaller the value the more likely to be missing
Missinggivenx2<-function(X){ 
  
  M<-rep(0, length(X))
  M <- rbinom(n=length(X), size=1, prob=1-exp(X)/(1+exp(X)))
  
  return(M)
  
}


M<-Missinggivenx2(X)


hist(X[M==0])

mean(X[M==0])



# 
# ###Adapted from Python code:
# 
# library(MASS)
# library(matrixStats)
# 
# # Dimension setting
# d <- 20
# 
# # Data generation function
# datagen <- function(n, theta, contamination = matrix(0, nrow = 0, ncol = d)) {
#   datag <- mvrnorm(n - nrow(contamination), mu = theta, Sigma = diag(d))
#   if (nrow(contamination) > 0) {
#     datag <- rbind(datag, contamination)
#   }
#   return(datag)
# }
# 
# # Create a matrix of differences between vector components
# outerdiff <- function(v1, v2) {
#   m <- v1 - t(v2)
#   return(m)
# }
# 
# # MMD estimation function
# MMD <- function(data, Nstep = 1000, gamma = 1, eta = 0.1) {
#   theta <- rep(0, d)
#   n <- nrow(data)
#   
#   for (t in 1:Nstep) {
#     for (j in 1:d) {
#       sample <- rnorm(n, theta[j], 1)
#       K1 <- outerdiff(sample, data[, j])
#       K2 <- outerdiff(sample, sample)
#       
#       diff <- theta[j] - sample
#       gradient <- -2 * (
#         mean(diff * (exp(-gamma * K2^2) - diag(n)) / (n - 1)) - 
#           mean(diff * exp(-gamma * K1^2))
#       )
#       
#       theta[j] <- theta[j] - eta * gradient / sqrt(t)
#     }
#   }
#   
#   return(theta)
# }
# 
# # Maximum Likelihood Estimation
# MLE <- function(data) {
#   return(colMeans(data))
# }
# 
# # Geometric Median
# GEO_MED <- function(data) {
#   return(colMedians(data))
# }
# 
# # Main simulation
# set.seed(123)  # For reproducibility
# theta <- rep(2, d)
# n <- 200
# N <- 20
# outl <- 9
# 
# resultats <- matrix(0, nrow = 3 * outl, ncol = N)
# 
# start_time <- Sys.time()
# 
# for (n_cont in 1:outl) {
#   for (i in 1:N) {
#     # Generate contamination
#     contamination <- matrix(rcauchy(5 * (n_cont - 1) * d), 
#                             ncol = d, 
#                             nrow = max(0, 5 * (n_cont - 1)))
#     
#     X <- datagen(n = n, theta = theta, contamination = contamination)
#     
#     resultats[3 * (n_cont - 1) + 1, i] <- mean((theta - MMD(X, gamma = d))^2)
#     resultats[3 * (n_cont - 1) + 2, i] <- mean((theta - MLE(X))^2)
#     resultats[3 * (n_cont - 1) + 3, i] <- mean((theta - GEO_MED(X))^2)
#   }
# }
# 
# vec_res <- sqrt(rowMeans(resultats))
# end_time <- Sys.time()
# 
# # Print results
# cat('Time:', as.numeric(end_time - start_time), 'seconds\n')
# 
# for (n_cont in 1:outl) {
#   cat('SqMSE of MMD for', (n_cont - 1) * 2.5, '% of outliers:', vec_res[3 * (n_cont - 1) + 1], '\n')
#   cat('SqMSE of MLE for', (n_cont - 1) * 2.5, '% of outliers:', vec_res[3 * (n_cont - 1) + 2], '\n')
#   cat('SqMSE of MED for', (n_cont - 1) * 2.5, '% of outliers:', vec_res[3 * (n_cont - 1) + 3], '\n')
# }