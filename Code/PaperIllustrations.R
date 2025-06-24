library(regMMD)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)



# Fonction pour calculer les quartiles de l'estimateur des extrêmes (Adversarial)
calculate_extreme_mean_quartiles_adversarial <- function(n, num_repetitions, alpha, epsilon) {
  extreme_estimators <- numeric(num_repetitions)
  
  for (i in 1:num_repetitions) {
    X <- rnorm(n, 0, 1)
    M <- rbinom(n, 1, alpha)
    
    
    num_smallest<-as.integer(epsilon/alpha*sum(M==0))-1
    #num_largest <- as.integer(n * epsilon / 2)
    #num_smallest <- as.integer(n * epsilon / 2)
    
    #indices_largest <- order(X, decreasing = TRUE)[1:num_largest]
    indices_smallest <- order(X)[1:num_smallest]
    
    #M[indices_largest] <- 0
    M[indices_smallest] <- 1
    
    
    
    
    # Estimation des extrêmes (maximum et minimum)
    max_value <- max(X[M == 0])
    min_value <- min(X[M == 0])
    
    # Moyenne des extrêmes
    extreme_estimators[i] <- (max_value + min_value) / 2
  }
  
  extreme_val <- mean(extreme_estimators)
  lower_quartile <- quantile(extreme_estimators, 0.25)
  upper_quartile <- quantile(extreme_estimators, 0.75)
  
  return(c(extreme_val, lower_quartile, upper_quartile))
}


# Fonction pour calculer les quartiles de l'estimateur des extrêmes (Adversarial)
calculate_extreme_mean_quartiles_realizable <- function(n, num_repetitions, alpha, epsilon) {
  extreme_estimators <- numeric(num_repetitions)
  
  for (i in 1:num_repetitions) {
    n1 <- as.integer((1 - epsilon) * n)
    X1 <- rnorm(n1, 0, 1)
    M1 <- rbinom(n1, 1, alpha)
    
    n2 <- n - n1
    X2 <- rnorm(n2, 0, 1)
    M2 <- 1-as.integer(X2 > 0)
    
    X <- c(X1, X2)
    M <- c(M1, M2)
    
    # Estimation des extrêmes (maximum et minimum)
    max_value <- max(X[M == 0])
    min_value <- min(X[M == 0])
    
    # Moyenne des extrêmes
    extreme_estimators[i] <- (max_value + min_value) / 2
  }
  
  extreme_val <- mean(extreme_estimators)
  lower_quartile <- quantile(extreme_estimators, 0.25)
  upper_quartile <- quantile(extreme_estimators, 0.75)
  
  return(c(extreme_val, lower_quartile, upper_quartile))
}



# Fonction pour calculer les quartiles de l'estimateur Mean (cas adversarial)
calculate_mean_quartiles_adversarial <- function(n, num_repetitions, alpha, epsilon) {
  mean_estimators <- numeric(num_repetitions)
  
  for (i in 1:num_repetitions) {
    X <- rnorm(n, 0, 1)
    M <- rbinom(n, 1, alpha)
    
    num_smallest<-as.integer(epsilon/alpha*sum(M==0))-1
    #num_largest <- as.integer(n * epsilon / 2)
    #num_smallest <- as.integer(n * epsilon / 2)
    
    #indices_largest <- order(X, decreasing = TRUE)[1:num_largest]
    indices_smallest <- order(X)[1:num_smallest]
    
    #M[indices_largest] <- 0
    M[indices_smallest] <- 1
    
    # Estimateur Mean
    mean_estimators[i] <- mean(X[M == 0])
  }
  
  mean_val <- mean(mean_estimators)
  lower_quartile <- quantile(mean_estimators, 0.25)
  upper_quartile <- quantile(mean_estimators, 0.75)
  
  return(c(mean_val, lower_quartile, upper_quartile))
}

# Fonction pour calculer les quartiles de l'estimateur Mean (cas réalisable)
calculate_mean_quartiles_realizable <- function(n, num_repetitions, alpha, epsilon) {
  mean_estimators <- numeric(num_repetitions)
  
  for (i in 1:num_repetitions) {
    n1 <- as.integer((1 - epsilon) * n)
    X1 <- rnorm(n1, 0, 1)
    M1 <- rbinom(n1, 1, alpha)
    
    n2 <- n - n1
    X2 <- rnorm(n2, 0, 1)
    M2 <- as.integer(X2 > 0)
    
    X <- c(X1, X2)
    M <- c(M1, M2)
    
    # Estimateur Mean
    mean_estimators[i] <- mean(X[M == 0])
  }
  
  mean_val <- mean(mean_estimators)
  lower_quartile <- quantile(mean_estimators, 0.25)
  upper_quartile <- quantile(mean_estimators, 0.75)
  
  return(c(mean_val, lower_quartile, upper_quartile))
}

# Fonction pour calculer les quartiles de l'estimateur MMD (cas réalisable)
calculate_mmd_quartiles_realizable <- function(n, num_repetitions, alpha, epsilon) {
  mmd_estimators <- numeric(num_repetitions)
  
  for (i in 1:num_repetitions) {
    n1 <- as.integer((1 - epsilon) * n)
    X1 <- rnorm(n1, 0, 1)  # Données réalisables
    M1 <- rbinom(n1, 1, alpha)
    
    n2 <- n - n1
    X2 <- rnorm(n2, 0, 1)  # Données restantes
    M2 <- 1-as.integer(X2 > 0)  # Indicateur réalisable
    
    X <- c(X1, X2)
    M <- c(M1, M2)
    
    # Estimateur MMD
    Est1 <- mmd_est(X[M == 0], model = "Gaussian.loc", par2 = 1)
    mmd_estimators[i] <- Est1$estimator
  }
  
  mean_val <- mean(mmd_estimators)
  lower_quartile <- quantile(mmd_estimators, 0.25)
  upper_quartile <- quantile(mmd_estimators, 0.75)
  
  return(c(mean_val, lower_quartile, upper_quartile))
}


# Fonction pour calculer les quartiles de l'estimateur MMD (cas adversarial)
calculate_mmd_quartiles_adversarial <- function(n, num_repetitions, alpha, epsilon) {
  mmd_estimators <- numeric(num_repetitions)
  
  for (i in 1:num_repetitions) {
    X <- rnorm(n, 0, 1)  # Génération des données adversariales
    M <- rbinom(n, 1, alpha)
    
    num_smallest<-as.integer(epsilon/alpha*sum(M==0))-1
    #num_largest <- as.integer(n * epsilon / 2)
    #num_smallest <- as.integer(n * epsilon / 2)
    
    #indices_largest <- order(X, decreasing = TRUE)[1:num_largest]
    indices_smallest <- order(X)[1:num_smallest]
    
    #M[indices_largest] <- 0
    M[indices_smallest] <- 1
    
    # Estimateur MMD
    Est1 <- mmd_est(X[M == 0], model = "Gaussian.loc", par2 = 1)
    mmd_estimators[i] <- Est1$estimator
  }
  
  mean_val <- mean(mmd_estimators)
  lower_quartile <- quantile(mmd_estimators, 0.25)
  upper_quartile <- quantile(mmd_estimators, 0.75)
  
  return(c(mean_val, lower_quartile, upper_quartile))
}




# Paramètres
sample_sizes <- seq(50, 1000, by = 50)  # Tailles d'échantillon
num_repetitions <- 500  # Nombre de répétitions pour chaque taille
alpha <- 0.5
epsilon <- 0.1

# Calcul des moyennes et des quartiles pour chaque taille d'échantillon
results_adversarial <- t(sapply(sample_sizes, calculate_mean_quartiles_adversarial, num_repetitions = num_repetitions, alpha = alpha, epsilon = epsilon))
results_realizable <- t(sapply(sample_sizes, calculate_mean_quartiles_realizable, num_repetitions = num_repetitions, alpha = alpha, epsilon = epsilon))

# Préparer les données pour le graphique
dfMean <- data.frame(
  SampleSize = rep(sample_sizes, 2),
  MeanEstimator = c(results_adversarial[, 1], results_realizable[, 1]),
  CI_Lower = c(results_adversarial[, 2], results_realizable[, 2]),
  CI_Upper = c(results_adversarial[, 3], results_realizable[, 3]),
  Setting = rep(c("Adversarial", "Realizable"), each = length(sample_sizes))
)



# Calcul des moyennes et des quartiles pour chaque taille d'échantillon
results_adversarial <- t(sapply(sample_sizes, calculate_mmd_quartiles_adversarial, num_repetitions = num_repetitions, alpha = alpha, epsilon = epsilon))
results_realizable <- t(sapply(sample_sizes, calculate_mmd_quartiles_realizable, num_repetitions = num_repetitions, alpha = alpha, epsilon = epsilon))

# Préparer les données pour le graphique
dfMMD <- data.frame(
  SampleSize = rep(sample_sizes, 2),
  MMD_Estimator = c(results_adversarial[, 1], results_realizable[, 1]),
  CI_Lower = c(results_adversarial[, 2], results_realizable[, 2]),
  CI_Upper = c(results_adversarial[, 3], results_realizable[, 3]),
  Setting = rep(c("Adversarial", "Realizable"), each = length(sample_sizes))
)





results_adversarial <- t(sapply(sample_sizes, calculate_extreme_mean_quartiles_adversarial, num_repetitions = num_repetitions, alpha = alpha, epsilon = epsilon))
results_realizable <- t(sapply(sample_sizes, calculate_extreme_mean_quartiles_realizable, num_repetitions = num_repetitions, alpha = alpha, epsilon = epsilon))

# Préparer les données pour le graphique
dfExtreme <- data.frame(
  SampleSize = rep(sample_sizes, 2),
  Extreme_Estimator = c(results_adversarial[, 1], results_realizable[, 1]),
  CI_Lower = c(results_adversarial[, 2], results_realizable[, 2]),
  CI_Upper = c(results_adversarial[, 3], results_realizable[, 3]),
  Setting = rep(c("Adversarial", "Realizable"), each = length(sample_sizes))
)

# Tracer les courbes avec zones de confiance basées sur les quartiles


##Plotting



# Create the three plots separately and assign them to variables
p1 <- ggplot(dfMean, aes(x = SampleSize, y = MeanEstimator, color = Setting, fill = Setting)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, linetype = 0) +
  labs(
    title = "MLE",
    x = "n",
    y = "",
    color = "Setting",
    fill = "Setting"
  ) +
  scale_y_continuous(
    limits = c(-0.25, 1),
    breaks = c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    #legend.position = "top",
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  )

p2 <- ggplot(dfMMD, aes(x = SampleSize, y = MMD_Estimator, color = Setting, fill = Setting)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, linetype = 0) +
  labs(
    title = "MMD",
    x = "n",
    y = "",
    color = "Setting",
    fill = "Setting"
  ) +
  scale_y_continuous(
    limits = c(-0.3, 1),
    breaks = c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)
  )+
  theme_minimal(base_size = 14) +
  theme(
    #legend.position = "top",
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  )

p3 <- ggplot(dfExtreme, aes(x = SampleSize, y = Extreme_Estimator, color = Setting, fill = Setting)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, linetype = 0) +
  labs(
    title = "Average of Extremes",
    x = "n",
    y = "",
    color = "Setting",
    fill = "Setting"
  ) +
  scale_y_continuous(
    limits = c(-0.25, 1),
    breaks = c(-0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9)
  )+
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    plot.margin = margin(20, 20, 20, 20)
  )

# Combine plots side-by-side using patchwork
combined_plot <- p1 + p2 + p3 + 
  plot_layout(nrow = 1, guides = "collect") +  # Collect legends
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  )

# Display the combined plot
combined_plot

# Save the combined plot as a PNG file
ggsave(
  filename = "combined_estimators.png", 
  plot = combined_plot, 
  width = 18, 
  height = 6, 
  dpi = 300
)








