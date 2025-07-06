##Gaussian MCAR with double contamination
---------------------------------------------
  
  
  
  # R Setup for Using https://github.com/cullena20/RobustMeanEstimation/tree/main
  # =======================================================

# Step 1: Install and load required R packages
# ---------------------------------------------
library(drf)
library(devtools)
library(reticulate)





# Step 2: Set up Python environment
# ----------------------------------

# Check if Python is available
py_available()

py_require(c("torch", "torchvision", "torchaudio"))
py_require(c("numpy", "scipy", "matplotlib", "seaborn", "pandas", "scikit-learn"))



# Step 5: Import Python modules and setup
# ----------------------------------------

# Import necessary Python modules
np <- import("numpy")
torch <- import("torch")
sys <- import("sys")
os <- import("os")
matplotlib <- import("matplotlib")
#eigenvalue_pruning<-import("eigenvalue_pruning")


source_python('RobustMeanEstimation-main/Algorithms/eigenvalue_pruning.py')
source_python('RobustMeanEstimation-main/Algorithms/que.py')
source_python('RobustMeanEstimation-main/Algorithms/que_utils.py')

############################################
# =======================================================




# Load required libraries
library(MASS)  # for multivariate normal generation
library(ggplot2)
library(reshape2)
source("mmd_est.R")
source("helpers.R")

# Set parameters
set.seed(123)  # for reproducibility
n <- 500      # total number of samples
d <- 10         # dimension (you can change this)
eps <- 0.3 # contamination fraction (10% contaminated data)


# Clean distribution parameters
mu_clean <- rep(0, d)        # mean vector of zeros
Sigma_clean <- diag(d)       # identity covariance matrix





Xtest<-createdistMCAR(n=5000, d=d, eps=eps, contdist="dirac", c=10)

sum(complete.cases(Xtest))/nrow(Xtest)

#mmd_estown(Xtest, model = "multidim.Gaussian.loc", par2 = 1, control=list(method="SGD", burnin=50, nstep=100))
##expected fraction:
## p2 / (1 - p1 + p2) + (p1 * (1 - p1)) / (d * (1 - p1 + p2)) * (1 - (p1 - p2)^d) / (1 - (p1 - p2))



contaminationlist<-c("rnorm", "rnorm", "rnorm", "dirac","dirac")
parameters<-c(0,1,10,1,10)
methodsp<-c("MMD", "mean", "median")
methodsnp<-c("evpruning", "que")
methods<-c(methodsp,methodsnp)

# Load required packages
library(foreach)
library(doParallel)

# Setup parallel backend for Windows
n_cores <- detectCores() - 1  # Use all cores except 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)


clusterExport(cl, c("contaminationlist", "parameters", 
                    "n", "d", "eps",  # Add your actual variable names here
                    "createdistMNARMCAR", "mmd_estown"))   # Add your function names here


clusterEvalQ(cl, {
  # library(your_package_name)  # Add any packages your functions need
  library(mice)
})

B<-50

#for (b in 1:B){
resultsp <- foreach(b = 1:B, .combine ='rbind') %dopar% {
  set.seed(b)
  
  cat(b)
  
  methodlist<-list()
  
  for (method in methodsp){
    methodlist[[method]] <- matrix(NA, nrow=1, ncol=length(contaminationlist))
    colnames(methodlist[[method]]) <- paste0(method, ".", contaminationlist,"_", parameters)
  }
  
  
  for (l in 1:length(contaminationlist)){
    
    contdist <- contaminationlist[l]
    
    cat(contdist)
    
    if (contdist=="rnorm"){
      X<-createdistMCAR(n=n, d=d, eps=eps, contdist=contdist, mu=parameters[l])
    } else if (contdist=="dirac"){
      
      X<-createdistMCAR(n=n, d=d, eps=eps, contdist=contdist, c=parameters[l]) 
      
    }
    
    for (method in methodsp){
      
      estimate<-get_estimate(X, method=method)
      methodlist[[method]][1,l]<-sqrt(sum(c(estimate)^2)) 
      
    }
  }
  
  res<-data.frame(do.call(cbind,methodlist))
  
  return(res)
  #return(list(MSEours=MSEours, MSEmean=MSEmean  ))
  
  
}



resultsnp <- t(sapply(1:B, function(b){
  set.seed(b)
  
  cat(b)
  
  methodlist<-list()
  
  for (method in methodsnp){
    methodlist[[method]] <- matrix(NA, nrow=1, ncol=length(contaminationlist))
    colnames(methodlist[[method]]) <- paste0(method, ".", contaminationlist,"_", parameters)
  }
  
  
  for (l in 1:length(contaminationlist)){
    
    contdist <- contaminationlist[l]
    
    cat(contdist)
    
    if (contdist=="rnorm"){
      X<-createdistMCAR(n=n, d=d, eps=eps, contdist=contdist, mu=parameters[l])
    } else if (contdist=="dirac"){
      
      X<-createdistMCAR(n=n, d=d, eps=eps, contdist=contdist, c=parameters[l]) 
      
    }
    
    
    for (method in methodsnp){
      
      estimate<-get_estimate(X, method=method)
      methodlist[[method]][1,l]<-sqrt(sum(c(estimate)^2)) 
      
    }
    
  }
  
  #res<-data.frame(cbind(MSEours, MSEmean, MSEmedian,MSEnormimp, MSEmeanimp ))
  
  res<-data.frame(do.call(cbind,methodlist))
  
  return(res)
  #return(list(MSEours=MSEours, MSEmean=MSEmean  ))
  
  
}))


###Need to concatenate resultsp and resultsnp!!!
resultsnp<-apply(resultsnp,2, unlist)
results<-cbind(resultsp, resultsnp)
####


save(results, file=paste0("MCAR_Results_", "eps_", eps,".R" ) )

resultsmean<-colMeans(results)
resultssd<-apply(results,2, sd)

methodlist<-list()
methodlistsd<-list()
for (method in methods){
  methodlist[[method]] <- rbind(resultsmean[grepl(paste0("^", method, "\\."), names(resultsmean))], resultssd[grepl(paste0("^", method, "\\."), names(resultssd))])
  #methodlistsd[[method]] <- resultssd[grepl(paste0("^", method, "\\."), names(resultssd))]
}




results_matrix <- round(data.frame(do.call(rbind, methodlist)),2)

# Create a copy for formatting
formatted_matrix <- results_matrix

# Find minimum value in each column and bold it
for(j in 1:ncol(results_matrix)) {
  # Add brackets to SD
  formatted_matrix[seq(from=2,to=(nrow(results_matrix)), by=2  ), j] <- paste0("(", results_matrix[seq(from=2,to=(nrow(results_matrix)), by=2  ), j], ")")
  min_idx <- which.min(results_matrix[seq(from=1,to=(nrow(results_matrix)-1), by=2  ), j])
  ##correction because we only look at odd numbers:
  min_idx <- min_idx + (min_idx-1) 
  formatted_matrix[min_idx, j] <- paste0("\\textbf{", results_matrix[min_idx, j], "}")
  formatted_matrix[min_idx+1, j] <- paste0("\\textbf{(", results_matrix[min_idx+1, j], ")}")
}

rownames(formatted_matrix)<-NULL
rownames(formatted_matrix)[seq(from=1,to=(nrow(results_matrix)-1), by=2  )] <-methods


# Convert to LaTeX table using xtable
library(xtable)
latex_table <- xtable(formatted_matrix, 
                      caption = "Comparison of MMD and MLE Results",
                      label = "tab:results")

# Print the LaTeX code
print(latex_table, 
      type = "latex", 
      include.rownames = TRUE,
      caption.placement = "top",
      sanitize.text.function = function(x) x)  # Don't escape LaTeX commands
