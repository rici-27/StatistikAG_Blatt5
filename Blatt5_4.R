## Aufgabe d)

n <- 100
get_ratio <- function(p, n=100, M=1000){
  error_sample_cov <- 0
  error_non_oracle <- 0
  Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
  Sigma_inv <- solve(Sigma)
  for (i in (1:M)){
    sample_cov <- get_sample_cov(n,p)
    non_oracle_est <- get_non_oracle_estimator(n,p)$estimator
    
    error_sample_cov <- error_sample_cov + (1/p) * norm(Sigma_inv - solve(sample_cov), type = "F")^2
    error_non_oracle <- error_non_oracle + (1/p) * norm(Sigma_inv - solve(non_oracle_est), type = "F")^2
  }
  ratio <- error_non_oracle/error_sample_cov
  return(ratio)
}

dim <- (1:10) * 10
ratios <- array(0, dim = 10)

for (i in (1:10)){
  ratios[i] <- get_ratio(p=dim[i])
  print(i)
}

View(ratios)

ggplot() + 
  geom_point(aes(x = dim, y = ratios)) 

############# hier nochmal schauen sieht bei daniel anders aus :_)

