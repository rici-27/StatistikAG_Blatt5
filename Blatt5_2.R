
# Aufgabe b)


w <- 0:100/100
gamma <- (1/5) * sum(diag(Sigma))
p <- 5

fehler_berechnen <- function(M, n, w){
  
  # Zwischenspeicher für Squared Bias, Variance, MSE
  sum_bias <- 0
  sum_var <- 0
  sum_mse <- 0
  
  # Einheitsmatrix für Shrinkage Schätzer
  I <- diag(1, nrow = p, ncol = p)
  
  for (i in (1:M)){
    X <- mvrnorm(n, mu, Sigma)
    S <- (1/n) * t(X) %*% X
    shrinkage <- w * gamma * I + (1-w) * S
    
    ############ hier weiter machen
    expectation_of_shrink_est <- (1-weight)* sigma_task_0 + weight*gamma*diag(dim)
    
    A <- expectation_of_shrink_est - sigma_task_0  #same for all j in M, not dependent on Sn!!
    #! might be calculated ouit of for loop
    sum_bias <- sum_bias + sum(diag(A %*% t(A)))
    
    B <- shrink_est - expectation_of_shrink_est #depends on j and Sn
    sum_var <- sum_var + sum(diag(B %*% t(B)))
    
    C <- shrink_est - sigma_task_0 #depends on j and Sn
    sum_mse <- sum_mse + sum(diag(C %*% t(C)))
    
    
  }
  return(storage)
}

  
# Data Frame Speicher für Fehlertherme und zugehöriges Gewicht
storage <- data.frame(array(0, dim = c(length(w), 4)))
colnames(storage) <- c("Weight", "Squared Bias", "Variance", "MSE")
View(storage)


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