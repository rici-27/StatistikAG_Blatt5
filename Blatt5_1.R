library(tidyverse)
library
library(MASS)

# Vorbereitung

n <- 100
Sigma <- matrix(c( 1 , rep (.1 ,4) , .1 , 1 , rep ( .1 , 3 ) , .1 , .1 , 1 , .1 , .1 , rep(.1 ,3) ,1 ,.1 , rep (.1 ,4) ,1) ,5 ,5)
mu <- c(rep(0,5))
View(Sigma)
?mvrnorm

X <- mvrnorm(n, mu, Sigma)
# X in R^{n \times p}
View(X)

## Aufgabe a)

M <- 1000
true_eigenvalues <- eigen(Sigma)$values
print(true_eigenvalues)

storage_ev <- matrix(0, nrow = M, ncol = 5)
for (i in (1:M)){
  X <- mvrnorm(n, mu, Sigma)
  S <- (1/n) * t(X) %*% X
  storage_ev[i,] <- eigen(S)$values
}
calc_ev <- sort((1/M) * colSums(storage_ev))
calc_ev

indizes <- c(1,2,3,4,5)

ggplot() +
  geom_point(aes(x=indizes, y = true_eigenvalues), size = 2, color = "blue") +
  geom_point(aes(x=indizes, y = calc_ev), size = 2, color = "darkgreen") +
  labs(x = "Indizes", y = "Eigenwerte", titel = "Wahre & Berechnete Eigenwerte")


