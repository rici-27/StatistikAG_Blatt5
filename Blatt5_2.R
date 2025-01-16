library(tidyverse)
library(MASS)

## Aufgabe b)

# Parameter
p <- 5
n <- 100
M <- 1000
Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
mu <- c(rep(0,p))
w <- 0:100/100
# Parameter gamma = (1/p) * trace(Sigma)
gamma <- (1/p) * sum(diag(Sigma))

# Funktion erstellen um Fehler-Terme für jeweiliges Gewicht w zu erhalten
# Parameter: Anzahl von Iterationen M, Beobachtungen n, Gewicht w
get_errors <- function(M, n, w){
  
  # Zwischenspeicher für Squared Bias, Variance, MSE
  sum_bias <- 0
  sum_var <- 0
  sum_mse <- 0
  
  for (i in (1:M)){
    
    X_transposed <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X_transposed) %*% X_transposed
    shrink_est <- w * gamma * diag(p) + (1-w) * S
    
    # Erwartungswert des Shrinkage Estimators (S ist erwartungstreu)
    expectation_of_shrink_est <- w * gamma * diag(p) + (1-w) * Sigma 
    
    # Squared Bias # gleich für jede Iteration!
    A <- expectation_of_shrink_est - Sigma
    sum_bias <- sum_bias + (1/p) * norm(A, type = "F")^2
    
    # Variance
    B <- shrink_est - expectation_of_shrink_est
    sum_var <- sum_var + (1/p) * norm(B, type = "F")^2
    
    # Mean Squared Error
    C <- shrink_est - Sigma
    sum_mse <- sum_mse + (1/p) * norm(C, type = "F")^2
    
  }
  
  squared_bias <- (1/M) * sum_bias
  variance <- (1/M) * sum_var
  mse <- (1/M) * sum_mse
  
  return(c(w, squared_bias, variance, mse))
}

  
# Data Frame als Speicher für Fehlerterme und zugehöriges Gewicht
storage <- data.frame(array(NA, dim = c(length(w), 4)))
colnames(storage) <- c("Weight", "SquaredBias", "Variance", "MSE")

# Berechnung
for (i in (1:length(w))){
  storage[i,] <- get_errors(M, n, w[i])
  print(i)
}

# Plot
ggplot() + 
  geom_point(aes(x=storage$Weight, y = storage$SquaredBias, color = "Squared Bias"),
             show.legend = TRUE) + 
  geom_point(aes(x=storage$Weight, y = storage$Variance, color = "Variance"),
             show.legend = TRUE) + 
  geom_point(aes(x=storage$Weight, y = storage$MSE, color = "MSE"),
             show.legend = TRUE) + 
  labs(x = "Gewicht w", y = "Werte",
       title = "Oracle Shrinkage Estimator Performance nach Gewicht, n = 100") +
  scale_color_manual(
    values = c("Squared Bias" = "darkblue", "Variance" = "darkgreen", "MSE" = "darkred"),
    name = "Legende"
  ) + 
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Empirisch optimales Gewicht w bestimmen
w_optimal <- storage$Weight[which.min(storage$MSE)]
w_optimal


# Durchführung mit n = 1000
storage_mod <- data.frame(array(0, dim = c(length(w), 4)))
colnames(storage_mod) <- c("Weight", "SquaredBias", "Variance", "MSE")

for (i in (1:length(w))){
  storage_mod[i,] <- get_errors(M, n=1000, w[i])
  print(i)
}

ggplot() + 
  geom_point(aes(x=storage_mod$Weight, y = storage_mod$SquaredBias,
                 color = "Squared Bias"), show.legend = TRUE) + 
  geom_point(aes(x=storage_mod$Weight, y = storage_mod$Variance,
                 color = "Variance"), show.legend = TRUE) + 
  geom_point(aes(x=storage_mod$Weight, y = storage_mod$MSE,
                 color = "MSE"), show.legend = TRUE) + 
  labs(x = "Gewicht w", y = "Werte",
       title = "Oracle Shrinkage Estimator Performance nach Gewicht, n = 1000") +
  scale_color_manual(
    values = c("Squared Bias" = "darkblue", "Variance" = "darkgreen", "MSE" = "darkred"),
    name = "Fehler"
  ) + 
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Empirisch optimales Gewicht w bestimmen
w_optimal_mod <- storage_mod$Weight[which.min(storage_mod$MSE)]
w_optimal_mod

