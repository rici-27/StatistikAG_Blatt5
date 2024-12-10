library(tidyverse)
library(MASS)

# Vorbereitung

n <- 100
Sigma <- matrix(c( 1 , rep (.1 ,4) , .1 , 1 , rep ( .1 , 3 ) , .1 , .1 , 1 , .1 , .1 , rep(.1 ,3) ,1 ,.1 , rep (.1 ,4) ,1) ,5 ,5)
mu <- c(rep(0,5))

## Aufgabe c)

M <- 1000
n <- 100

# Pilot Estimator = Sample Covariance Matrix
# soll für das schätzen vom optimal weight auch schon monte carlo gemacht werden?
# falls ja diese funktion verwenden, aber welches X?
# erstmal lassen

get_pilot <- function(M, n){
  storage <- array(0, dim = c(1000, 5, 5))
  I <- diag(1, nrow = 5, ncol = 5)
  for (i in (1:M)){
    X <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X) %*% X
    storage[i,,] <- S
  }
  mean_scm <- apply(storage, c(2,3), mean)

  return(mean_scm)
}

# Pilot Estimator Teil 1)

X <- mvrnorm(n, mu, Sigma)
sample_cov <- (1/n) * t(X) %*% X


# Parameter schätzen

# gamma
gamma <- sum(diag(sample_cov)) 

# delta^2
I <- diag(1, nrow = 5, ncol = 5)
delta2 <- norm((1/5) * (sample_cov - gamma * I), type = "F")^2

# beta^2
summe <- sum(sapply(1:n, function(k) {
  X_k <- X[k,]                     
  diff <- X_k %*% t(X_k) - sample_cov                
  norm(diff, type = "F")^2                
}))
beta2_pilot <- (1/n^2) * summe
beta2 <- min(beta2_pilot, delta2)

# alpha^2
alpha2 <- delta2 - beta2

# jetzt Gewichte berechnen
rho1 <- (beta2/delta2) * gamma
rho2 <- alpha2/delta2

get_estimator <- function(M, w1, w2){
  storage <- array(0, dim = c(1000, 5, 5))
  I <- diag(1, nrow = 5, ncol = 5)
  for (i in (1:M)){
    X <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X) %*% X
    estimator <- w1 * I + w2 * S
    storage[i,,] <- estimator
  }
  mean_estimator <- apply(storage, c(2,3), mean)
  return(mean_estimator)
}

non_oracle <- get_estimator(M,rho1, rho2)

# Eigenwerte vergleichen

true_eigenvalues <- eigen(Sigma)$values
estimated_ev <- eigen(non_oracle)$values

ggplot() +
  geom_point(aes(x=indizes, y = true_eigenvalues), size = 2, color = "blue") +
  geom_point(aes(x=indizes, y = estimated_ev), size = 2, color = "darkgreen") +
  labs(x = "Indizes", y = "Eigenwerte", titel = "Wahre & Berechnete Eigenwerte")

# Optimale Gewichte anschauen
rho1
rho2

# rho2 ist nah dran, und rho1 = w * gamma circa, also passt ungefähr
# komisch dass ergebnisse dann nicht besser sind



