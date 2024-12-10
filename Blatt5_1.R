library(tidyverse)
library(MASS)

# Vorbereitung

n <- 100
Sigma <- matrix(c( 1 , rep (.1 ,4) , .1 , 1 , rep ( .1 , 3 ) , .1 , .1 , 1 , .1 , .1 , rep(.1 ,3) ,1 ,.1 , rep (.1 ,4) ,1) ,5 ,5)
mu <- c(rep(0,5))
#View(Sigma)
?mvrnorm

X <- mvrnorm(n, mu, Sigma)
# X in R^{n \times p}
#View(X)

## Aufgabe a)

M <- 1000
true_eigenvalues <- sort(eigen(Sigma)$values)
# , symmetric=TRUE, only.values=TRUE vll noch benutzen
print(true_eigenvalues)

storage_ev <- matrix(0, nrow = M, ncol = 5)
for (i in (1:M)){
  X <- mvrnorm(n, mu, Sigma)
  S <- (1/n) * t(X) %*% X
  storage_ev[i,] <- sort(eigen(S)$values)
  # sort in welche richtung?
}
calc_ev <- sort((1/M) * colSums(storage_ev))
calc_ev

?sort
indizes <- c(1,2,3,4,5)

ggplot() +
  geom_point(aes(x=indizes, y = true_eigenvalues), size = 2, color = "blue") +
  geom_point(aes(x=indizes, y = calc_ev), size = 2, color = "darkgreen") +
  labs(x = "Indizes", y = "Eigenwerte", titel = "Wahre & Berechnete Eigenwerte") +
  theme_minimal()

calc_ev_df <- data.frame(Eigenwerte = c(storage_ev))

ggplot(calc_ev_df, aes(x=Eigenwerte)) +
  geom_histogram(fill="pink", binwidth = 0.05) +
  geom_vline(xintercept = mean(true_eigenvalues), color = "blue", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = mean(c(calc_ev)), color = "red", linetype = "dotted", linewidth = 1) +
  theme_minimal()
# legen und überschrift fehlt usw.

# vllt auch sum in histogram packen

# Aufgabe b)


w <- 0:100/100
gamma <- (1/5) * sum(diag(Sigma))
p <- 5

fehler_berechnen <- function(M, n, w){
  storage <- array(0, dim = c(1000, 5, 5))
  I <- diag(1, nrow = 5, ncol = 5)
  for (i in (1:M)){
    X <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X) %*% X
    shrinkage <- w * gamma * I + (1-w) * S
    storage[i,,] <- shrinkage
  }
  mean_shrinkage <- apply(storage, c(2,3), mean)
  # hier nochmal frobenius norm prüfen
  squared_bias <- (1/p)*(norm(mean_shrinkage - Sigma, type="F"))^2
  variance <- (1/p)*norm(mean_shrinkage - (w * gamma * I + (1-w) * Sigma), type = "F")^2
  mse <- (1/p)*norm(mean_shrinkage - Sigma, type = "F")^2
  return(c(squared_bias, variance, mse))
}

ergebnisse <- data.frame(matrix(nrow = 101, ncol = 4))
colnames(ergebnisse) <- c('w', 'SquaredBias', 'Variance', 'MSE')
ergebnisse$w <- w


for (i in (1:101)){
  ergebnisse[i,2:4] <- fehler_berechnen(M, n=100, w[i])
  print(i)
}
View(ergebnisse)

# hier stimmt was gar nicht

ggplot() + 
  geom_point(aes(x=ergebnisse$w, y = ergebnisse$SquaredBias, color = "Squared Bias"), show.legend = TRUE) + 
  geom_point(aes(x=ergebnisse$w, y = ergebnisse$Variance, color = "Variance"), show.legend = TRUE) + 
  geom_point(aes(x=ergebnisse$w, y = ergebnisse$MSE, color = "MSE"), show.legend = TRUE) + 
  labs(x = "Gewicht w", y = "Werte", title = "Shrinkage Estimator Performance je nach Gewicht, n = 100") +
  scale_color_manual(
    values = c("Squared Bias" = "darkblue", "Variance" = "darkgreen", "MSE" = "darkred"),
    name = "Fehler"
  ) + 
  #ylim(0, 0.0005) +
  theme_classic()


# empirical optimal w ? kleinster mse?
index <- which.min(ergebnisse$MSE)
w_optimal <- ergebnisse$w[index]
w_optimal

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

index <- which.min(ergebnisse_mod$MSE)
w_optimal_mod <- ergebnisse$w[index]
w_optimal_mod

# optimales gewicht 0 ist bissi wild, passt der Erwartungswert?