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
  storage_ev[i,] <- sort(eigen(S)$values)
}
calc_ev <- sort((1/M) * colSums(storage_ev))
calc_ev

indizes <- c(1,2,3,4,5)

ggplot() +
  geom_point(aes(x=indizes, y = true_eigenvalues), size = 2, color = "blue") +
  geom_point(aes(x=indizes, y = calc_ev), size = 2, color = "darkgreen") +
  labs(x = "Indizes", y = "Eigenwerte", titel = "Wahre & Berechnete Eigenwerte")

calc_ev_df <- data.frame( Eigenwerte = c(storage_ev))

ggplot(calc_ev_df, aes(x=Eigenwerte)) +
  geom_histogram(fill="pink", binwidth = 0.05) +
  geom_vline(xintercept = mean(true_eigenvalues), color = "blue", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = mean(c(calc_ev)), color = "red", linetype = "dotted", linewidth = 1) 
  


# Aufgabe b)


w <- 0:100/100
gamma <- (1/5) * sum(diag(Sigma))

fehler_berechnen <- function(M, n, w){
  storage <- matrix(0, nrow = 1000, ncol = 3)
  I <- diag(1, nrow = 5, ncol = 5)
  for (i in (1:M)){
    X <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X) %*% X
    shrinkage <- w * gamma * I + (1-w) * S
    # bias
    storage[i,1] <- (1/5) * norm(shrinkage - Sigma, type = "F")
    # variance 
    storage[i,2] <- ((1/5) * norm(shrinkage - (w * gamma * I + (1-w) * Sigma), type = "F"))^2
    # MSE
    storage[i,3] <- ((1/5) * norm(shrinkage - Sigma))^2
  }
  squared_bias <- mean(storage[i,1])^2
  variance <- mean(storage[i,2])
  mse <- mean(storage[i,3])
  return(c(squared_bias, variance, mse))
}

ergebnisse <- data.frame(matrix(nrow = 101, ncol = 4))
colnames(ergebnisse) <- c('w', 'SquaredBias', 'Variance', 'MSE')
ergebnisse$w <- w

for (i in (1:101)){
  ergebnisse[i,2:4] <- fehler_berechnen(M, n=100, w[i])
}
View(ergebnisse)

ggplot() + 
  geom_point(aes(x=ergebnisse$w, y = ergebnisse$SquaredBias, color = "Squared Bias"), show.legend = TRUE) + 
  geom_point(aes(x=ergebnisse$w, y = ergebnisse$Variance, color = "Variance"), show.legend = TRUE) + 
  geom_point(aes(x=ergebnisse$w, y = ergebnisse$MSE, color = "MSE"), show.legend = TRUE) + 
  labs(x = "Gewicht w", y = "Werte", title = "Shrinkage Estimator Performance je nach Gewicht, n = 100") +
  scale_color_manual(
    values = c("Squared Bias" = "darkblue", "Variance" = "darkgreen", "MSE" = "darkred"),
    name = "Fehler"
  ) + 
  theme_classic()

# jetzt mit n = 1000

ergebnisse_mod <- data.frame(matrix(nrow = 101, ncol = 4))
colnames(ergebnisse_mod) <- c('w', 'SquaredBias', 'Variance', 'MSE')
ergebnisse_mod$w <- w

for (i in (1:101)){
  ergebnisse_mod[i,2:4] <- fehler_berechnen(M, n=1000, w[i])
  print(i)
}

ggplot() + 
  geom_point(aes(x=ergebnisse_mod$w, y = ergebnisse_mod$SquaredBias, color = "Squared Bias"), show.legend = TRUE) + 
  geom_point(aes(x=ergebnisse_mod$w, y = ergebnisse_mod$Variance, color = "Variance"), show.legend = TRUE) + 
  geom_point(aes(x=ergebnisse_mod$w, y = ergebnisse_mod$MSE, color = "MSE"), show.legend = TRUE) + 
  labs(x = "Gewicht w", y = "Werte", title = "Shrinkage Estimator Performance je nach Gewicht, n = 1000") +
  scale_color_manual(
    values = c("Squared Bias" = "darkblue", "Variance" = "darkgreen", "MSE" = "darkred"),
    name = "Fehler"
  ) + 
  theme_classic()