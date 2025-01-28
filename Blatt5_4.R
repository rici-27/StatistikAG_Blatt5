library(tidyverse)
library(MASS)

## Aufgabe d)

n <- 100
M <- 1000
get_ratio <- function(p, n=100, M=1000){
  error_sample_cov <- 0
  error_non_oracle <- 0
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
  Sigma_inv <- solve(Sigma)
  
  for (i in (1:M)){
    estimators_list <- get_estimators(n, p)
    sample_cov <- estimators_list$sample_cov
    non_oracle_est <- estimators_list$shrinkage
    
    error_sample_cov <- error_sample_cov +
                        (1/p) * norm(Sigma_inv - solve(sample_cov), type = "F")^2
    error_non_oracle <- error_non_oracle +
                        (1/p) * norm(Sigma_inv - solve(non_oracle_est), type = "F")^2
  }
  
  # Division durch M für mittleren Fehler
  error_non_oracle <- error_non_oracle/M
  error_sample_cov <- error_sample_cov/M
  ratio <- error_non_oracle/error_sample_cov
  return(list(ratio=ratio,
              error_sample_cov = error_sample_cov,
              error_non_oracle = error_non_oracle))
}

dim <- (1:10) * 10
# Array als Speicher für Ratios
ratios <- array(NA, dim = length(dim))
# Array als Speicher für Fehler
# Spalte 1 für Sample Cov, Spalte 2 für Non Oracle Estimator
errors <- array(NA, dim = c(length(dim), 2))


for (i in (1:length(dim))){
  ratio_list <- get_ratio(p=dim[i])
  ratios[i] <- ratio_list$ratio
  errors[i,1] <- ratio_list$error_sample_cov
  errors[i,2] <- ratio_list$error_non_oracle
  print(i)
}

# Plot für die Ratios
ggplot() + 
  geom_point(aes(x = dim, y = ratios)) +
  labs(x = "Dimension", y = "Ratio",
       title = "Vergleich der Schätzfehler von Non Oracle und Sample Cov") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Plot für die Fehler
ggplot() + 
  geom_point(aes(x= dim[1:5], y = errors[(1:5),1], color = "Sample Cov"), size = 2) +
  geom_point(aes(x= dim[1:5], y = errors[(1:5),2], color = "Non Oracle"), size = 2) +
  labs(x = "Dimension", y = "Fehler",
       title = "Vergleich der Schätzfehler von Non Oracle und Sample Cov") +
  scale_color_manual(
    name = "Legende",
    values = c("Sample Cov" = "purple", "Non Oracle" = "orange")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )
