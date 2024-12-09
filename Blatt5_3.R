## Aufgabe d)

n <- 100
get_sample_cov <- function(M, p){
  storage <- array(0, dim = c(1000, p, p))
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, p, p)
  mu <- c(rep(0,p))
  for (i in (1:M)){
    X <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X) %*% X
    storage[i,,] <- S
  }
  mean_sample_cov <- apply(storage, c(2,3), mean)
  return(mean_sample_cov)
}


get_shrinkage_est <- function(p){
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
}

get_ratio <- function(p){
  
}