

# X<-as.matrix(cbind(df$log_income_observed, df$age))
# 
# res<-mmd_estown(X, model = "multidim.Gaussian.loc", par2 = 1, control=list(method="SGD"))
# 
# 
# ### TO DO for package:
# ## Simply reweigh the dist function output u in the kernel functions K1d, Kmd
# ## Ignore the test that checks that X is numeric
# ## Add par = median(x, na.rm=T)
# 











mmd_estown<-function (x, model, par1 = NULL, par2 = NULL, kernel = "Gaussian", 
          bdwth = "median", control = list()) {
            list.models.multidim = c("multidim.Gaussian.loc", "multidim.Gaussian.scale", 
                                     "multidim.Gaussian", "multidim.Dirac")
            list.models.uni = c("Gaussian.loc", "Gaussian.scale", "Gaussian", 
                                "continuous.uniform.loc", "continuous.uniform.upper", 
                                "continuous.uniform.lower.upper", "exponential", "gamma.shape", 
                                "gamma.rate", "gamma", "Cauchy", "Pareto", "beta", "Poisson", 
                                "geometric", "Dirac", "binomial.prob", "binomial.size", 
                                "binomial", "discrete.uniform")
            list.kernels = c("Gaussian", "Laplace", "Cauchy")
            list.methods = c("GD", "SGD", "EXACT")
            if (is.null(control$burnin) == FALSE) {
              if ((is.double(control$burnin) == FALSE) || (length(control$burnin) != 
                                                           1)) {
                res$error = c(res$error, "control$burnin must be numerical")
              }
              else if (floor(control$burnin) != control$burnin) {
                res$error = c(res$error, "control$burnin must be an integer")
              }
              else burnin = control$burnin
            }
            else burnin = 500
            if (is.null(control$nstep) == FALSE) {
              if ((is.double(control$nstep) == FALSE) || (length(control$nstep) != 
                                                          1)) {
                res$error = c(res$error, "control$nstep must be numerical")
              }
              else if (floor(control$nstep) != control$nstep) {
                res$error = c(res$error, "control$nstep must be an integer")
              }
              else nstep = control$nstep
            }
            else nstep = 1000
            if (is.null(control$epsilon) == FALSE) {
              if ((is.double(control$epsilon) == FALSE) || (length(control$epsilon) != 
                                                            1)) {
                res$error = c(res$error, "control$epsilon must be numerical")
              }
              else epsilon = control$epsilon
            }
            else epsilon = 10^(-5)
            if (is.null(control$stepsize) == TRUE) 
              stepsize = "auto"
            else stepsize = control$stepsize
            if (is.null(control$method) == TRUE) 
              method = "default"
            else method = control$method
            out = list(model = model, par1init = par1, par2init = par2, 
                       kernel = kernel, bdwth = bdwth, burnin = burnin, nstep = nstep, 
                       stepsize = stepsize, epsilon = epsilon, method = method, 
                       error = NULL, estimator = NULL, type = "est")
            if (model %in% list.models.uni) {
              #if ((is.vector(x) == FALSE) || (is.numeric(x) == FALSE)) 
              #  out$error = c(out$error, "Data x should be a numeric vector for this model")
            }
            else if (model %in% list.models.multidim) {
              if (is.matrix(x) == FALSE) 
                out$error = c(out$error, "Data x should be an (n,p) matrix for this model (n: sample size)")
            }
            else out$error = c(out$error, "Model unknown")
            if (is.null(out$error) == FALSE) {
              for (i in 1:length(out$error)) warning(out$error[i])
              return(invisible(out))
            }
            if (bdwth == "median") {
              ##Weirdly weighting the missing values:
              ## Original code:
              #bdwth = median(c(dist(x, method = "euclidean")))
              ### What we do in the algorithm:
              u = as.matrix(dist(x, method = "euclidean"))
              weighting<-sqrt(as.matrix(dist(is.na(rbind(x))*1,method = "manhattan"))/ncol(x))
              weighting[weighting==0] <- 1
              u<-(u*weighting)
              bdwth=median(c(u))
              ######
              if (bdwth == 0) 
                bdwth = 1
            }
            else if (is.double(bdwth) == FALSE) 
              out$error = c(out$error, "bdwth must be numeric or 'median'")
            else if (bdwth <= 0) 
              out$error = c(out$error, "bdwth must be positive")
            if (kernel %in% list.kernels == FALSE) 
              out$error = c(out$error, "Kernel unknown")
            if (is.double(burnin) == FALSE) 
              out$error = c(out$error, "burnin must be numeric")
            if (is.double(nstep) == FALSE) 
              out$error = c(out$error, "nstep must be numeric")
            if (is.double(epsilon) == FALSE) 
              out$error = c(out$error, "epsilon must be numeric")
            if (is.double(stepsize) == FALSE) 
              if (stepsize != "auto") 
                out$error = c(out$error, "stepsize must be numeric or 'auto'")
            if (is.null(out$error) == FALSE) {
              for (i in 1:length(out$error)) warning(out$error[i])
              return(invisible(out))
            }
            if (model == "Gaussian.loc") {
              if (((kernel != "Gaussian") && (method == "GD")) || (method == 
                                                                   "EXACT")) 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              if (((kernel == "Gaussian") && (method == "GD")) || ((kernel == 
                                                                    "Gaussian") && (method == "default"))) {
                method = "GD"
                resultat = GD.MMD.Gaussian.loc(x = x, par1 = par1, 
                                               par2 = par2, kernel = kernel, bdwth = bdwth, 
                                               burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                               epsilon = epsilon)
              }
              if ((method == "SGD") || ((kernel != "Gaussian") && (method == 
                                                                   "default"))) {
                method = "GD"
                resultat = SGD.MMD.Gaussian.locown(x = x, par1 = par1, 
                                                par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                epsilon = epsilon)
              }
            }
            if (model == "continuous.uniform.loc") {
              if (method == "EXACT") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              if ((method == "GD") || (method == "default")) {
                method = "GD"
                resultat = GD.MMD.continuous.uniform.loc(x = x, par1 = par1, 
                                                         par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                         burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                         epsilon = epsilon)
              }
              if (method == "SGD") 
                resultat = SGD.MMD.continuous.uniform.loc(x = x, 
                                                          par1 = par1, par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                          burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                          epsilon = epsilon)
            }
            if (model == "binomial.prob") {
              if (method == "EXACT") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              if ((method == "GD") || (method == "default")) {
                method = "GD"
                resultat = GD.MMD.binomial.prob(x = x, par1 = par1, 
                                                par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                epsilon = epsilon)
              }
              if (method == "SGD") 
                resultat = SGD.MMD.binomial.prob(x = x, par1 = par1, 
                                                 par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                 burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                 epsilon = epsilon)
            }
            if (model == "multidim.Gaussian.loc") {
              if (((kernel != "Gaussian") && (method == "GD")) || (method == 
                                                                   "EXACT")) 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              if (((kernel == "Gaussian") && (method == "GD")) || ((kernel == 
                                                                    "Gaussian") & (method == "default"))) {
                method = "GD"
                resultat = GD.MMD.multidim.Gaussian.loc(x = x, par1 = par1, 
                                                        par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                        burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                        epsilon = epsilon)
              }
              if ((method == "SGD") || ((kernel != "Gaussian") && (method == 
                                                                   "default"))) {
                method = "SGD"
                resultat = SGD.MMD.multidim.Gaussian.locown(x = x, par1 = par1, 
                                                         par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                         burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                         epsilon = epsilon)
              }
            }
            if (model == "Gaussian.scale") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.Gaussian.scale(x = x, par1 = par1, 
                                                  par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                  burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                  epsilon = epsilon)
              }
            }
            if (model == "Gaussian") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.Gaussian(x = x, par1 = par1, par2 = par2, 
                                            kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                            nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "multidim.Gaussian") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.multidim.Gaussian(x = x, par1 = par1, 
                                                     par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                     burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                     epsilon = epsilon)
              }
            }
            if (model == "multidim.Gaussian.scale") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.multidim.Gaussian.scale(x = x, 
                                                           par1 = par1, par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                           burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                           epsilon = epsilon)
              }
            }
            if (model == "continuous.uniform.upper") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.continuous.uniform.upper(x = x, 
                                                            par1 = par1, par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                            burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                            epsilon = epsilon)
              }
            }
            if (model == "continuous.uniform.lower.upper") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.continuous.uniform.lower.upper(x = x, 
                                                                  par1 = par1, par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                                  burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                                  epsilon = epsilon)
              }
            }
            if (model == "exponential") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.gamma.rate(x = x, par1 = 1, par2 = par1, 
                                              kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                              nstep = nstep, stepsize = stepsize, epsilon = epsilon)
                resultat$par1 = resultat$par2
                resultat$par2 = NULL
              }
            }
            if (model == "gamma.shape") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.gamma.shape(x = x, par1 = par1, 
                                               par2 = par2, kernel = kernel, bdwth = bdwth, 
                                               burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                               epsilon = epsilon)
              }
            }
            if (model == "gamma.rate") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.gamma.rate(x = x, par1 = par1, 
                                              par2 = par2, kernel = kernel, bdwth = bdwth, 
                                              burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                              epsilon = epsilon)
              }
            }
            if (model == "gamma") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.gamma(x = x, par1 = par1, par2 = par2, 
                                         kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                         nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "Cauchy") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.Cauchy(x = x, par1 = par1, par2 = par2, 
                                          kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                          nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "Pareto") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.Pareto(x = x, par1 = par1, par2 = par2, 
                                          kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                          nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "beta") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.beta(x = x, par1 = par1, par2 = par2, 
                                        kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                        nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "Poisson") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.Poisson(x = x, par1 = par1, par2 = par2, 
                                           kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                           nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "geometric") {
              if (method == "default") 
                method = "SGD"
              if (method != "SGD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = SGD.MMD.geometric(x = x, par1 = par1, 
                                             par2 = par2, kernel = kernel, bdwth = bdwth, 
                                             burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                             epsilon = epsilon)
              }
            }
            if (model == "Dirac") {
              if (method == "default") 
                method = "GD"
              if (method != "GD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = GD.MMD.Dirac(x = x, par1 = par1, par2 = par2, 
                                        kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                        nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "multidim.Dirac") {
              if (method == "default") 
                method = "GD"
              if (method != "GD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = GD.MMD.multidim.Dirac(x = x, par1 = par1, 
                                                 par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                 burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                 epsilon = epsilon)
              }
            }
            if (model == "binomial") {
              if (method == "default") 
                method = "GD"
              if (method != "GD") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = GD.MMD.binomial(x = x, par1 = par1, par2 = par2, 
                                           kernel = kernel, bdwth = bdwth, burnin = burnin, 
                                           nstep = nstep, stepsize = stepsize, epsilon = epsilon)
              }
            }
            if (model == "binomial.size") {
              if (method == "default") 
                method = "EXACT"
              if (method != "EXACT") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = EXACT.MMD.binomial.size(x = x, par1 = par1, 
                                                   par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                   burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                   epsilon = epsilon)
              }
            }
            if (model == "discrete.uniform") {
              if (method == "default") 
                method = "EXACT"
              if (method != "EXACT") 
                out$error = c(out$error, "This method is not possible with this kernel and model")
              else {
                resultat = EXACT.MMD.discrete.uniform(x = x, par1 = par1, 
                                                      par2 = par2, kernel = kernel, bdwth = bdwth, 
                                                      burnin = burnin, nstep = nstep, stepsize = stepsize, 
                                                      epsilon = epsilon)
              }
            }
            if (is.null(out$error) == FALSE) {
              for (i in 1:length(out$error)) warning(out$error[i])
              return(invisible(out))
            }
            if (is.null(resultat$error) == FALSE) {
              out$error = resultat$error
              for (i in 1:length(out$error)) warning(out$error[i])
              return(invisible(out))
            }
            else {
              out$par1init = resultat$par1
              out$par2init = resultat$par2
              out$stepsize = resultat$stepsize
              out$bdwth = resultat$bdwth
              out$error = resultat$error
              out$method = method
              out$estimator = resultat$estimator
              class(out) <- "estMMD"
              return(invisible(out))
            }
          }



K1down <- function (x, y, kernel = "Laplace", bdwth = 1) 
{
  x[is.na(y)]<-0
  y[is.na(y)]<-0
  
  x[is.na(x)]<-0
  y[is.na(x)]<-0
  
  u = outer(x, y, FUN = "-")/bdwth
  if (kernel == "Gaussian") {
    return(exp(-u^2))
  }
  else if (kernel == "Laplace") {
    return(exp(-abs(u)))
  }
  else if (kernel == "Cauchy") {
    return(1/(2 + u^2))
  }
}


Kmdown<-function (x, y, kernel = "Laplace", bdwth = 1) 
{
  

  n = dim(x)[1]
  
  ###Careful the dist function is weighting the observations with missing values!!
  ## https://stackoverflow.com/questions/18117174/function-dist-not-behaving-as-expected-on-vectors-with-missing-values
  u = as.matrix(dist(rbind(x, y), diag = TRUE, upper = TRUE, 
                     method = "euclidean"))[1:n, (n + 1):(2 * n)]

  ##With missing values we need to correct the reweighting done in the dist function
  weighting<-sqrt(as.matrix(dist(is.na(rbind(x, y))*1, diag = TRUE, upper = TRUE,method = "manhattan"))/ncol(x))[1:n, (n + 1):(2 * n)]
  weighting[weighting==0] <- 1
  
  u<-(u*weighting)/bdwth
  
  
  if (kernel == "Gaussian") {
    return(exp(-u^2))
  }
  else if (kernel == "Laplace") {
    return(exp(-abs(u)))
  }
  else if (kernel == "Cauchy") {
    return(1/(2 + u^2))
  }
}



SGD.MMD.Gaussian.locown<-function (x, par1, par2, kernel, bdwth, burnin, nstep, stepsize, 
          epsilon) 
{
  n = length(x)
  res = list(par1 = par1, par2 = par2, stepsize = stepsize, 
             bdwth = bdwth, error = NULL, estimator = NULL)
  if (is.null(par1)) {
    par = median(x, na.rm=T)
  }
  else if ((is.double(par1)) && (length(par1) == 1)) {
    par = par1
  }
  else {
    res$error = c(res$error, "par1 must be numerical")
  }
  if (is.null(par2)) {
    res$error = c(res$error, "par2 missing")
  }
  else if ((is.double(par2) == FALSE) || (length(par2) != 1)) {
    res$error = c(res$error, "par2 must be numerical")
  }
  else if (par2 <= 0) {
    res$error = c(res$error, "par2 must be positive")
  }
  if (is.null(res$error) == FALSE) 
    return(res)
  if (stepsize == "auto") 
    stepsize = par2
  norm.grad = epsilon
  res$par1 = par
  res$par2 = par2
  res$stepsize = stepsize
  for (i in 1:burnin) {
    x.sampled = rnorm(n = n, mean = par, sd = par2)
    ker = (K1down(x.sampled, x.sampled, kernel = kernel, bdwth = bdwth) - 
             diag(n))/(n - 1) - K1down(x.sampled, x, kernel = kernel, 
                                    bdwth = bdwth)/n
    gradL = (x.sampled - par)/(par2^2)
    grad = 2 * mean(gradL %*% ker)
    norm.grad = norm.grad + grad^2
    par = par - stepsize * grad/sqrt(norm.grad)
  }
  par_mean = par
  for (i in 1:nstep) {
    x.sampled = rnorm(n = n, mean = par, sd = par2)
    ker = (K1down(x.sampled, x.sampled, kernel = kernel, bdwth = bdwth) - 
             diag(n))/(n - 1) - K1down(x.sampled, x, kernel = kernel, 
                                    bdwth = bdwth)/n
    gradL = (x.sampled - par)/(par2^2)
    grad = 2 * mean(gradL %*% ker)
    norm.grad = norm.grad + grad^2
    par = par - stepsize * grad/sqrt(norm.grad)
    par_mean = (par_mean * i + par)/(i + 1)
  }
  res$estimator = par_mean
  return(res)
}


SGD.MMD.multidim.Gaussian.locown<-function (x, par1, par2, kernel, bdwth, burnin, nstep, stepsize, 
          epsilon) 
{
  n = dim(x)[1]
  p = dim(x)[2]
  res = list(par1 = par1, par2 = par2, stepsize = stepsize, 
             bdwth = bdwth, error = NULL, estimator = NULL)
  if (is.null(par1)) {
    par = rep(0, p)
    for (i in 1:p) par[i] = median(x[, i], na.rm=T)
  }
  else if ((is.vector(par1) == FALSE) || (is.numeric(par1) == 
                                          FALSE)) {
    res$error = c(res$error, "par1 must be a numerical vector")
  }
  else if (length(par1) != p) {
    res$error = c(res$error, "wrong dimension for par1")
  }
  else {
    par = par1
  }
  if (is.null(par2)) {
    res$error = c(res$error, "par2 missing")
  }
  else if ((is.double(par2) == FALSE) || (length(par2) != 1)) {
    res$error = c(res$error, "par2 must be numerical")
  }
  else if (par2 <= 0) {
    res$error = c(res$error, "par2 must be positive")
  }
  if (is.null(res$error) == FALSE) 
    return(res)
  if (stepsize == "auto") 
    stepsize = par2
  norm.grad = epsilon
  res$par1 = par
  res$par2 = par2
  res$stepsize = stepsize
  
  
  tmp<-order_by_missingness_pattern(x)
  
  x<-tmp$ordered_data
  M<- is.na(x) * 1
  patterns<-unique(M)
  print(paste0("Number of Patterns: ", nrow(patterns)))
  
  for (i in 1:burnin) {
    #cat(i)
    x.sampled.centered = matrix(data = rnorm(n = n * p, mean = 0, 
                                             sd = par2), nrow = n, ncol = p)
    x.sampled = x.sampled.centered + matrix(data = par, nrow = n, 
                                            ncol = p, byrow = TRUE)
    ##Need to add NAs, This is not correct!!!!!
    #x.sampled.NA<-x.sampled
    #x.sampled.NA[is.na(x)]<-NA
    
    kerXY<-Kmdown(x.sampled, x, kernel = kernel, 
                  bdwth = bdwth)/n
    

    gradL = x.sampled.centered/(par2^2)
    
    
    
  #  gradYY<-matrix(NA, nrow=nrow(patterns), ncol=p)
  # for (i in 1:nrow(patterns)) {
  # 
  #     pattern <- patterns[i,]
  #     # Find rows with this pattern
  #     pattern_rows <- which(rowSums(unname(M==matrix(pattern, nrow = nrow(M), ncol = p, byrow = TRUE)))==p)
  # 
  #     np<-length(pattern_rows)
  # 
  #     x.sampled.NA<-x.sampled
  #     x.sampled.NA[matrix(is.na(x[pattern_rows,][1,]),nrow=n, ncol=p, byrow = T)]<-NA
  # 
  #     kerXYi = (Kmdown(x.sampled, x.sampled.NA, kernel = kernel, bdwth = bdwth) -
  #                diag(n))/(n - 1)
  # 
  # 
  #     gradYY[i,]<-c(2 * t(gradL) %*% kerXYi %*% matrix(data = 1/n,
  #                                         nrow = n, ncol = 1))
  # 
  # }
    
    gradYY<-t(apply(patterns, 1, function(pattern) {      
      pattern_rows <- which(rowSums(unname(M==matrix(pattern, nrow = nrow(M), ncol = p, byrow = TRUE)))==p)
    
    np<-length(pattern_rows)
    
    x.sampled.NA<-x.sampled
    x.sampled.NA[matrix(is.na(x[pattern_rows,][1,]),nrow=n, ncol=p, byrow = T)]<-NA
    
    kerXYi = (Kmdown(x.sampled, x.sampled.NA, kernel = kernel, bdwth = bdwth) - 
                diag(n))/(n - 1) 
    
    
   c(2 * t(gradL) %*% kerXYi %*% matrix(data = 1/n, 
                                                     nrow = n, ncol = 1))} ))
    

    gradYY= colSums(gradYY)/n
    gradXY = c(2 * t(gradL) %*% kerXY %*% matrix(data = 1/n, 
                                             nrow = n, ncol = 1))
    
    grad<-gradYY - gradXY
    
    
    norm.grad = norm.grad + sum(grad^2)
    par = par - stepsize * grad/sqrt(norm.grad)
  }
  par_mean = par
  for (i in 1:nstep) {
    #cat(i)
    x.sampled.centered = matrix(data = rnorm(n = n * p, mean = 0, 
                                             sd = par2), nrow = n, ncol = p)
    x.sampled = x.sampled.centered + matrix(data = par, nrow = n, 
                                            ncol = p, byrow = TRUE)
    
    ##Need to add NAs (This step is not clear yet!)
    gradL = x.sampled.centered/(par2^2)
    
    # gradYY<-matrix(NA, nrow=nrow(patterns), ncol=p)
    # for (i in 1:nrow(patterns)) {
    #   
    #   pattern <- patterns[i,]
    #   # Find rows with this pattern
    #   pattern_rows <- which(rowSums(unname(M==matrix(pattern, nrow = nrow(M), ncol = p, byrow = TRUE)))==p)
    #   
    #   np<-length(pattern_rows)
    #   
    #   x.sampled.NA<-x.sampled
    #   x.sampled.NA[matrix(is.na(x[pattern_rows,][1,]),nrow=n, ncol=p, byrow = T)]<-NA
    #   
    #   kerXYi = (Kmdown(x.sampled, x.sampled.NA, kernel = kernel, bdwth = bdwth) - 
    #               diag(n))/(n - 1) 
    #   
    #   
    #   gradYY[i,]<-c(2 * t(gradL) %*% kerXYi %*% matrix(data = 1/n, 
    #                                                    nrow = n, ncol = 1))
    #   
    # }
    
    
    gradYY<-t(apply(patterns,1, function(pattern){
      pattern_rows <- which(rowSums(unname(M==matrix(pattern, nrow = nrow(M), ncol = p, byrow = TRUE)))==p)
      
      np<-length(pattern_rows)
      
      x.sampled.NA<-x.sampled
      x.sampled.NA[matrix(is.na(x[pattern_rows,][1,]),nrow=n, ncol=p, byrow = T)]<-NA
      
      kerXYi = (Kmdown(x.sampled, x.sampled.NA, kernel = kernel, bdwth = bdwth) - 
                  diag(n))/(n - 1) 
      
      
     c(2 * t(gradL) %*% kerXYi %*% matrix(data = 1/n, 
                                                       nrow = n, ncol = 1))}  ))
    
    gradYY= colSums(gradYY)/n
    
    
    kerXY<-Kmdown(x.sampled, x, kernel = kernel, 
                  bdwth = bdwth)/n
    gradXY = c(2 * t(gradL) %*% kerXY %*% matrix(data = 1/n, 
                                                 nrow = n, ncol = 1))
    
    grad<-gradYY - gradXY
    
    
      norm.grad = norm.grad + sum(grad^2)
    par = par - stepsize * grad/sqrt(norm.grad)
    par_mean = (par_mean * i + par)/(i + 1)
  }
  res$estimator = par_mean
  return(res)
}


order_by_missingness_pattern <- function(X) {
  # Convert matrix to data frame if it's not already
  if (is.matrix(X)) {
    X_df <- as.data.frame(X)
  } else {
    X_df <- X
  }
  
  # Create a binary matrix where 1 = missing, 0 = not missing
  miss_matrix <- is.na(X_df) * 1
  
  # Create a pattern ID for each row
  # This converts the binary pattern to a single number that uniquely identifies each pattern
  pattern_id <- numeric(nrow(miss_matrix))
  for (i in 1:ncol(miss_matrix)) {
    pattern_id <- pattern_id + miss_matrix[, i] * 2^(i-1)
  }
  
  # Order rows by pattern ID (lowest first = fewest missing values)
  ordered_idx <- order(pattern_id)
  
  # Return the ordered data
  return(list(
    ordered_data = X_df[ordered_idx, ],
    ordered_idx = ordered_idx,
    patterns = pattern_id[ordered_idx]
  ))
}

