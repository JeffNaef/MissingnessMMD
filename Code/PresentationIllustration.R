

source("mmd_est.R")


# Function to generate data with missing values for income and age
generate_data_with_missing_values <- function(n_samples = 1000, 
                                              mu_income = 10.5, 
                                              mu_age = 40, 
                                              sigma_income = 1.0, 
                                              sigma_age = 15, 
                                              correlation = 0.3,
                                              alpha = 0.2, 
                                              epsilon = 0.3, 
                                              seed = 42,
                                              mech="MAR") {
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Define the mean vector and covariance matrix
  mu <- c(mu_income, mu_age)
  
  # Compute covariance from correlation and standard deviations
  cov_income_age <- correlation * sigma_income * sigma_age
  sigma <- matrix(c(
    sigma_income^2, cov_income_age,
    cov_income_age, sigma_age^2
  ), nrow = 2)
  
  # Generate samples from the multivariate normal distribution
  library(MASS)
  X <- mvrnorm(n = n_samples, mu = mu, Sigma = sigma)
  
  # Split into log_income and age
  log_income <- X[, 1]
  age <- X[, 2]
  
  # Generate missing pattern based on age
  # P[X₁ missing | X = x] = πₘ₂(x) = (1 - ε)α + ε1{x₂ > 50}
  if (mech=="MAR"){
  missing_prob <- (1 - epsilon) * alpha + epsilon * (age > 50)
  }else{
    missing_prob <- (1 - epsilon) * alpha + epsilon * (log_income > 11) 
  }
  # Sample missing indicators
  income_missing <- rbinom(n_samples, 1, missing_prob)
  income_observed <- as.logical(!income_missing)
  
  # Create patterns: m₁ = (0, 0) and m₂ = (1, 0)
  # m₁: both observed, m₂: only age observed
  pattern <- rep(0, n_samples)
  pattern[income_missing == 1] <- 2  # m₂ = (1, 0)
  pattern[income_missing == 0] <- 1  # m₁ = (0, 0)
  
  # Create masked version of log_income
  log_income_observed <- log_income
  log_income_observed[income_missing == 1] <- NA
  
  # Create actual income in dollars for interpretability
  income <- exp(log_income)
  income_observed <- exp(log_income_observed)
  
  # Create a data frame
  df <- data.frame(
    log_income = log_income,  # True log income (would be unknown in real dataset)
    log_income_observed = log_income_observed,  # Observed log income with NAs
    age = age,
    income_observed = income_observed,
    pattern = pattern,
    missing_prob = missing_prob,
    income = income,
    income_observed_dollars = income_observed
  )
  
  return(df)
}

# Function to visualize the missing data patterns
visualize_missing_data <- function(df) {
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  
  # Create a factor for observed/missing status
  df$missing_status <- factor(ifelse(is.na(df$log_income_observed), 
                                     "Income Missing", 
                                     "Income Observed"))
  
  # Plot 1: Scatter plot of age vs log income, colored by missing status
  p1 <- ggplot(df, aes(x = age, y = log_income, color = missing_status)) +
    geom_point(alpha = 0.5) +
    labs(title = "Age vs. Log Income (True Values)",
         x = "Age", 
         y = "Log Income") +
    theme_minimal() +
    #theme(legend.title = element_blank())
    theme(legend.position = "none")
    
  # Plot 2: Histogram of age by missing status
  p2 <- ggplot(df, aes(x = age, fill = missing_status)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 20) +
    labs(title = "Distribution of Age by Missing Status",
         x = "Age", 
         y = "Count") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Plot 1: Scatter plot of age vs log income, colored by missing status
  p3 <- ggplot(df, aes(x = age, y = log_income)) +
    geom_point(alpha = 0.5) +
    labs(title = "Age vs. Log Income (True Values)",
         x = "Age", 
         y = "Log Income") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Plot 4: Missing percentage by age group
  df$age_group <- cut(df$age, 
                      breaks = c(0, 30, 40, 50, 60, 70, 100),
                      labels = c("<30", "30-40", "40-50", "50-60", "60-70", ">70"))
  
  p4 <- df %>%
    group_by(age_group) %>%
    summarise(prop_missing = mean(is.na(log_income_observed))) %>%
    ggplot(aes(x = age_group, y = prop_missing)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Proportion of Missing Income by Age Group",
         x = "Age Group", 
         y = "Proportion Missing") +
    theme_minimal()
  
  # Arrange all plots in a grid
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}




epsilon<-0.2
alpha<-0.3

# Example usage
# Generate data with reasonable parameters for income and age
df <- generate_data_with_missing_values(
  n_samples = 2000,
  mu_income = 10.5,  # corresponds to about $36,000
  mu_age = 40,
  sigma_income = 0.8,
  sigma_age = 15,
  correlation = 0.3,  # slight positive correlation between age and income
  alpha = alpha,  # bas missing rate
  epsilon = epsilon,  # additional missing for older people
  mech="MNAR"
)

# Display summary statistics
mean(df[df$age > 50,]$pattern==2)
mean(df[df$age <= 50,]$pattern==2)

df[df$pattern==0,]

# Count of missing vs. observed
missing_count <- table(is.na(df$log_income_observed))
cat("\nMissing pattern counts:\n")
cat("Observed income:", missing_count["FALSE"], "\n")
cat("Missing income:", missing_count["TRUE"], "\n")
cat("Missing percentage:", missing_count["TRUE"] / sum(missing_count) * 100, "%\n")

# Visualize the data
visualize_missing_data(df)


df <- generate_data_with_missing_values(
  n_samples = 200,
  mu_income = 10.5,  # corresponds to about $36,000
  mu_age = 40,
  sigma_income = 0.8,
  sigma_age = 15,
  correlation = 0.3,  # slight positive correlation between age and income
  alpha = alpha,  # bas missing rate
  epsilon = epsilon,  # additional missing for older people
  mech="MAR"
)


X<-as.matrix(cbind(df$log_income_observed, df$age))


#X<-as.matrix(cbind(df$log_income, df$age))

res<-mmd_estown(X, model = "multidim.Gaussian.loc", par2 = 1, control=list(method="SGD"))

cat("MMD Estimator for MAR", "(", res$estimator, ")")
cat("Truth", "(", c(10.5,40), ")")

mean(df$log_income_observed,na.rm=T)



# Example usage
# Generate data with reasonable parameters for income and age
df <- generate_data_with_missing_values(
  n_samples = 200,
  mu_income = 10.5,  # corresponds to about $36,000
  mu_age = 40,
  sigma_income = 0.8,
  sigma_age = 15,
  correlation = 0.3,  # slight positive correlation between age and income
  alpha = alpha,  # base 20% missing rate
  epsilon = epsilon,  # additional 30% missing for older people,
  mech="MNAR"
)


X<-as.matrix(cbind(df$log_income_observed, df$age))


#X<-as.matrix(cbind(df$log_income, df$age))

res<-mmd_estown(X, model = "multidim.Gaussian.loc", par2 = 1, control=list(method="SGD"))

cat("MMD Estimator for MNAR", "(", res$estimator, ")")
cat("Truth", "(", c(10.5,40), ")")

mean(df$log_income_observed,na.rm=T)



######################
###Simulation study###
######################
library(foreach)
library(doParallel)

cores <- detectCores() - 1  # Use all cores except one
cl <- makeCluster(cores)
registerDoParallel(cl)

# Simulation study
B <- 500
mechanism<-"MNAR"

# Export necessary functions and variables to the cluster
clusterExport(cl, c("generate_data_with_missing_values", "mmd_estown", "alpha", "epsilon"))

# Parallelize the loop
MMDest <- foreach(b = 1:B, .combine = rbind) %dopar% {
  
  df <- generate_data_with_missing_values(
    n_samples = 200,
    mu_income = 10.5,  # corresponds to about $36,000
    mu_age = 40,
    sigma_income = 0.8,
    sigma_age = 15,
    correlation = 0.3,  # slight positive correlation between age and income
    alpha = alpha,  # base missing rate
    epsilon = epsilon,  # additional missing for older people
    mech = mechanism,
    seed = b
  )
  
  X <- as.matrix(cbind(df$log_income_observed, df$age))
  
  res <- mmd_estown(X, model = "multidim.Gaussian.loc", par2 = 1, control = list(method = "SGD"))
  res$estimator
}

# Stop the cluster when done
stopCluster(cl)




# After running your parallelized simulation
# MMDest now contains B rows with 2 columns (one for each estimator)

# Load ggplot2 for better visualization
library(ggplot2)

# Convert the matrix to a data frame for easier plotting
results_df <- as.data.frame(MMDest)
colnames(results_df) <- c("Income", "Age")

# True means
true_means <- c(10.5, 40)
est_mean<-colMeans(results_df)

# Create histogram for Income estimator
p1 <- ggplot(results_df, aes(x = Income)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "white", alpha = 0.7) +
  geom_vline(xintercept = est_mean[1], color = "black", linetype = "dashed", size = 1, alpha = 0.7) +
  geom_vline(xintercept = true_means[1], color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Income Estimator",
       subtitle = paste("True mean =", true_means[1], "(black)" , "Estimate =", round(est_mean[1],2), "(red)" ),
       x = "Estimated Value",
       y = "Frequency") +
  theme_minimal()



# Create histogram for Age estimator
p2 <- ggplot(results_df, aes(x = Age)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "white", alpha = 0.7) +
  geom_vline(xintercept = est_mean[2], color = "black", linetype = "dashed", size = 1, alpha = 0.7) +
  geom_vline(xintercept = true_means[2], color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Age Estimator",
       subtitle = paste("True mean =", true_means[2], "(black)" , "Estimate =", round(est_mean[2],2), "(red)" ),
       x = "Estimated Value", 
       y = "Frequency") +
  theme_minimal()

# Display plots side by side
library(gridExtra)
grid.arrange(p1, p2, ncol = 2)

combined_plot <- grid.arrange(p1, p2, ncol = 2)
if (mechanism=="MAR"){
ggsave("combined_estimator_histograms_MAR.png", combined_plot, width = 12, height = 5, dpi = 300)
}else{
  ggsave("combined_estimator_histograms_MNAR.png", combined_plot, width = 12, height = 5, dpi = 300)
}

