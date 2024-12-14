library(tidyverse)
library(MASS)

## Aufgabe c)

# Parameter
p <- 5
n <- 100
M <- 1000
Sigma <- matrix(c( 1 , rep (.1 ,4) , .1 , 1 , rep ( .1 , 3 ) , .1 , .1 , 1 , .1 , .1 , rep(.1 ,3) ,1 ,.1 , rep (.1 ,4) ,1) ,5 ,5)
mu <- c(rep(0,p))
w <- 0:100/100
gamma <- (1/5) * sum(diag(Sigma))


# Funktion zur Berechnung des Non-Oracle Estimator

get_non_oracle_estimator <- function(n, p){
  # Daten simulieren
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
  mu = c(rep(0,p))
  X <- mvrnorm(n, mu, Sigma)
  sample_cov <- (1/n) * t(X) %*% X
  
  # Parameter schätzen
  # gamma
  gamma <- (1/p) * sum(diag(sample_cov %*% diag(p)))

  # delta^2
  delta2 <- (1/5) * norm((sample_cov - gamma * diag(p)), type = "F")^2
  
  # beta^2
  beta2_help <- rep(0, n)
  for (i in 1:n){
    temp <- outer(t(X)[,i] , X[i,]) - sample_cov
    # Produkt zweier Vektoren um (n x n)-Matrix zu erhalten
    beta2_help[i] <- (1/p) * norm(temp, type = "F")^2 
  }
  beta2_pilot <- (1/n^2) * sum(beta2_help)
  beta2 <- min(beta2_pilot, delta2)
  
  # alpha^2
  alpha2 <- delta2 - beta2
  
  # Gewichte berechnen
  rho1 <- (beta2/delta2) * gamma
  rho2 <- alpha2/delta2
  
  # Schätzer berechen
  estimator <- rho1 * diag(p) + rho2 * sample_cov
  return(list(estimator = estimator, para1 = rho1, para2 = rho2))
}

get_sample_cov <- function(n, p){
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
  mu = c(rep(0,p))
  X <- mvrnorm(n, mu, Sigma)
  sample_cov <- (1/n) * t(X) %*% X
  return(sample_cov)
}


# Monte-Carlo Simulation und Berechnung der Eigenwerte

get_estimated_ev <- function(M, n, p){
  
  sum_ev_sample_cov <- rep(0,p)
  for (i in (1:M)){
    estimator <- get_sample_cov(n, p)
    sum_ev_sample_cov <- sum_ev_sample_cov + sort(eigen(estimator, symmetric=TRUE, only.values=TRUE)$values, decreasing = FALSE) 
  }
  
  sum_ev_non_oracle <- rep(0,p)
  for (i in (1:M)){
    estimator <- get_non_oracle_estimator(n, p)$estimator
    sum_ev_non_oracle <- sum_ev_non_oracle + sort(eigen(estimator, symmetric=TRUE, only.values=TRUE)$values, decreasing = FALSE) 
  }
  
  estimated_ev <- data.frame(matrix(0, p, 2))
  colnames(estimated_ev) <- c("SampleCov", "NonOracle")
  
  estimated_ev$SampleCov <- sum_ev_sample_cov/M
  estimated_ev$NonOracle <- sum_ev_non_oracle/M
  
  return(estimated_ev)
}

estimated_ev <- get_estimated_ev(1000, 100, 5)
true_eigenvalues <- sort(eigen(Sigma)$values)
indizes <- (1:5)

ggplot() +
  geom_point(aes(x = indizes, y = estimated_ev$SampleCov, color = "SampleCov", shape = "SampleCov"), size = 3) +
  geom_point(aes(x = indizes, y = estimated_ev$NonOracle, color = "NonOracle", shape = "NonOracle"), size = 3) +
  geom_point(aes(x = indizes, y = true_eigenvalues, color = "Sigma", shape = "Sigma"), size = 3) +
  labs(x = "Indizes", y = "Eigenwerte", title = "Vergleich wahrer und geschätzer Eigenwerte") + 
  scale_color_manual(
    values = c("SampleCov" = "darkblue", "NonOracle" = "darkgreen", "Sigma" = "darkred"),
    name = "Schätzer"
  ) + 
  scale_shape_manual(
    values = c("SampleCov" = 16, "NonOracle" = 17, "Sigma" = 18), # Form der Punkte
    name = "Schätzer"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )


# Alle Eigenwerte und Weights ausgeben lassen
get_all_est_ev_and_weights <- function(M, n, p){
  # 3-dim Array zur Speicherung aller Eigenwerte
  # Dimension 1 für Sample Cov, Dimension 2 für Non Oracle Estimator
  ev_array <- array(0, dim=c(2, M*p, 2))
  
  # Speicher für weights
  weights <- array(0, dim=c(1000,2))

  for (i in (1:M)){
    print(i)
    non_oracle <- get_non_oracle_estimator(n, p)
    non_oracle_est <- non_oracle$estimator
    weights[i,] <- c(non_oracle$para1, non_oracle$para2)
    sample_cov <- get_sample_cov(n, p)
    non_oracle_ev <- sort(eigen(non_oracle_est, symmetric=TRUE, only.values=TRUE)$values, decreasing = FALSE) 
    sample_cov_ev <- sort(eigen(sample_cov, symmetric=TRUE, only.values=TRUE)$values, decreasing = FALSE) 
    
    
    
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
colnames(all_ev_sample_cov) <- c("Index", "Wert")
all_ev_sample_cov$Schätzer <- "SampleCov"

all_ev_non_oracle <- data.frame(all_est_ev[2,,])
colnames(all_ev_non_oracle) <- c("Index", "Wert")
all_ev_non_oracle$Schätzer <- "NonOracle"

all_ev_df <- rbind(all_ev_non_oracle, all_ev_sample_cov)
all_ev_df$Index <- as.factor(all_ev_df$Index)

ggplot(all_ev_df, aes(x=Schätzer, y=Wert)) +
  geom_boxplot(aes(fill=Index)) +
  labs(title = "Vergleich von Non Oracle Estimator und Sample Cov") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

######### hier noch irgendwie weights plotten 

