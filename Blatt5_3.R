library(tidyverse)
library(MASS)

## Aufgabe c)

# Parameter
p <- 5
n <- 100
M <- 1000
Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
mu <- c(rep(0,p))
w <- 0:100/100


# Funktion zur Berechnung des Non-Oracle Estimator & Sample Covariance

get_estimators <- function(n, p){
  # Daten simulieren
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
  mu <- c(rep(0,p))
  X_transposed <- mvrnorm(n, mu, Sigma)
  sample_cov <- (1/n) * t(X_transposed) %*% X_transposed
  
  # Parameter schätzen
  # gamma
  gamma <- (1/p) * sum(diag(sample_cov %*% diag(p)))

  # delta^2
  delta2 <- (1/p) * norm((sample_cov - gamma * diag(p)), type = "F")^2
  
  # beta^2
  beta2_help <- rep(0, n)
  for (i in 1:n){
    temp <- outer(t(X_transposed)[,i] , X_transposed[i,]) - sample_cov
    # Produkt zweier Vektoren um (p x p)-Matrix zu erhalten
    beta2_help[i] <- (1/p) * norm(temp, type = "F")^2 
  }
  beta2_pilot <- (1/n^2) * sum(beta2_help)
  beta2 <- min(beta2_pilot, delta2)
  
  # alpha^2
  alpha2 <- delta2 - beta2
  
  # Gewichte berechnen
  para1 <- (beta2/delta2)
  para2 <- alpha2/delta2
  
  # Schätzer berechen
  shrinkage <- para1 * gamma * diag(p) + para2 * sample_cov
  return(list(shrinkage = shrinkage,
              sample_cov = sample_cov,
              para1 = para1,
              para2 = para2))
}


# Alle Eigenwerte und Weights ausgeben lassen
get_all_est_ev_and_weights <- function(M, n, p){
  # 3-dim Array zur Speicherung aller Eigenwerte
  # Dimension 1 für Sample Cov, Dimension 2 für Non Oracle Estimator
  ev_array <- array(NA, dim=c(2, M*p, 2))
  
  # Speicher für weights
  weights <- array(NA, dim=c(M, 2))

  for (i in (1:M)){
    estimators_list <- get_estimators(n, p)
    non_oracle_est <- estimators_list$shrinkage
    weights[i,] <- c(estimators_list$para1, estimators_list$para2)
    sample_cov <- estimators_list$sample_cov
    non_oracle_ev <- sort(eigen(non_oracle_est, symmetric=TRUE,
                                only.values=TRUE)$values, decreasing = FALSE) 
    sample_cov_ev <- sort(eigen(sample_cov, symmetric=TRUE,
                                only.values=TRUE)$values, decreasing = FALSE) 

    
    for (j in (1:p)){
      ev_array[1, (i-1)*p + j, 1] <- j
      ev_array[1, (i-1)*p + j, 2] <- sample_cov_ev[j]
      
      ev_array[2, (i-1)*p + j, 1] <- j
      ev_array[2, (i-1)*p + j, 2] <- non_oracle_ev[j]
    }
  }
  return(list(eigenvalues = ev_array, weights = weights))
}

all_iterations <- get_all_est_ev_and_weights(M=1000, n=100, p=5)
all_est_ev <- all_iterations$eigenvalues
weights <- all_iterations$weights

# Data Frames für BoxPlots erstellen

all_ev_sample_cov <- data.frame(all_est_ev[1,,])
colnames(all_ev_sample_cov) <- c("Index", "Eigenwert")
all_ev_sample_cov$Schätzer <- "SampleCov"

all_ev_non_oracle <- data.frame(all_est_ev[2,,])
colnames(all_ev_non_oracle) <- c("Index", "Eigenwert")
all_ev_non_oracle$Schätzer <- "NonOracle"

all_ev_df <- rbind(all_ev_non_oracle, all_ev_sample_cov)
all_ev_df$Index <- as.factor(all_ev_df$Index)

ggplot(all_ev_df, aes(x=Index, y=Eigenwert)) +
  geom_boxplot(aes(fill=Schätzer), outlier.colour = "black", outlier.size = 1) +
  geom_segment(aes(x = 0.5, xend = 4.5, y = 0.9, yend = 0.9, color = "Eigenwert 0.9"), 
               linetype = "dashed", size = 0.8, inherit.aes = FALSE) +
  geom_segment(aes(x = 4.5, xend = 5.5, y = 1.4, yend = 1.4, color = "Eigenwert 1.4"), 
               linetype = "dashed", size = 0.8, inherit.aes = FALSE) +
  scale_color_manual(
    values = c("Eigenwert 0.9" = "black", "Eigenwert 1.4" = "brown"),
    name = "Legende"
  ) +
  labs(title = "Vergleich der sortierten Eigenwerte von Non Oracle Estimator und Sample Cov") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Plot Highlight
ggplot(all_ev_df, aes(x=Schätzer, y=Eigenwert)) +
  geom_boxplot(aes(fill=Index)) +
  labs(title = "Vergleich der sortierten Eigenwerte von Non Oracle Estimator und Sample Cov") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Vergleich der geschätzen optimalen Gewichte zu den Ergebnissen aus b)
weights_df <- data.frame(weights[,1])
colnames(weights_df) = c("Optimal_Weight")
weights_df$Boxplot <- as.factor(" ")

ggplot(weights_df, aes(x = Boxplot, y = Optimal_Weight)) + 
  geom_boxplot(color = "purple", fill = "purple", alpha = 0.2,
               outlier.colour = "red", outlier.fill = "red", outlier.size = 2)  +
  geom_hline(aes(yintercept = w_optimal, color = "Gewicht für n=100"),
             linetype = "dashed") +
  geom_hline(aes(yintercept = w_optimal_mod, color = "Gewicht für n=1000"),
             linetype = "dashed") +
  scale_color_manual(
    name = "Referenzlinien", 
    values = c("Gewicht für n=100" = "blue", "Gewicht für n=1000" = "orange")
  ) +
  labs(title = "Vergleich der geschätzen optimalen Gewichte") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )
