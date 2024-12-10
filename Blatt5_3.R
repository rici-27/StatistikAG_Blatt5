## Aufgabe d)

n <- 100
get_sample_cov <- function(p, M){
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


get_shrinkage_est <- function(p, M){
  
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, p, p)
  mu <- c(rep(0,p))
  X <- mvrnorm(n, mu, Sigma)
  sample_cov <- (1/n) * t(X) %*% X
  
  # Parameter schätzen
  gamma <- sum(diag(sample_cov)) 
  I <- diag(1, nrow = p, ncol = p)
  delta2 <- norm((1/p) * (sample_cov - gamma * I), type = "F")^2
  summe <- sum(sapply(1:n, function(k) {
    X_k <- X[k,]                     
    diff <- X_k %*% t(X_k) - sample_cov                
    norm(diff, type = "F")^2                
  }))
  beta2_pilot <- (1/n^2) * summe
  alpha2 <- delta2 - beta2
  # jetzt Gewichte berechnen
  rho1 <- (beta2/delta2) * gamma
  rho2 <- alpha2/delta2

  storage <- array(0, dim = c(1000, p, p))
  
  for (i in (1:M)){
    X <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X) %*% X
    estimator <- rho1 * I + rho2 * S
    storage[i,,] <- estimator
  }
  mean_estimator <- apply(storage, c(2,3), mean)
  return(mean_estimator)
}

#get_sample_cov(p=5, M=1000)

#get_shrinkage_est(p=5, M=1000)


get_ratio <- function(p){
  print(p)
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, p, p)
  Sigma_inv <- solve(Sigma)
  fehler_sample_cov <- ((1/p) * norm(Sigma_inv- solve(get_sample_cov(p=p, M=1000)), type = "2"))^2
  fehler_shrinkage <- ((1/p) * norm(Sigma_inv - solve(get_shrinkage_est(p=p, M=1000)), type= "2"))^2
  ratio <- fehler_sample_cov / fehler_shrinkage
}


max_l <- 30
ratios <- array(0, dim = max_l)

for (p in (1:max_l)){
  ratios[p] <- get_ratio(p)
}

View(ratios)
p_werte <- (1:max_l)

ggplot() + 
  geom_point(aes(x = p_werte, y = ratios))


