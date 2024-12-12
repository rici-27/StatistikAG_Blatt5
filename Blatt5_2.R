
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