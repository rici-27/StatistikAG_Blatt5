library(tidyverse)
library(MASS)

# Vorbereitung

n <- 100
p <- 5
# Kovarianz Matrix und Mittelwert erstellen
Sigma <- matrix(c( 1 , rep (.1 ,4) , .1 , 1 , rep ( .1 , 3 ) , .1 , .1 , 1 , .1 , .1 , rep(.1 ,3) ,1 ,.1 , rep (.1 ,4) ,1) ,5 ,5)
mu <- c(rep(0,p))
# Multivariate Daten erzeugen
X <- mvrnorm(n, mu, Sigma)
# Vorsicht: X in R^{n \times p}


## Aufgabe a)

M <- 1000
# symmetric -> Matrix ist symmetrisch
# only.values -> Eigenvektoren werden nicht berechnet
# decreasing = FALSE -> Eigenwerte werden aufsteigend sortiert
true_eigenvalues <- sort(eigen(Sigma, symmetric=TRUE, only.values=TRUE)$values, decreasing = FALSE)

# Speicher für Eigenwerte erstellen
storage_ev <- matrix(0, nrow = M, ncol = p)
# Monte Carlo Simulation
for (i in (1:M)){
  X <- mvrnorm(n, mu, Sigma)
  S <- (1/n) * t(X) %*% X
  storage_ev[i,] <- sort(eigen(S)$values)
}
# Mittelwerte der Eigenwerte berechnen
calc_ev <- (1/M) * colSums(storage_ev)


# Boxplot
# erst Data Frame erstellen
eigenvalues_df <- data.frame(matrix(0, nrow=p*M, ncol=2))
colnames(eigenvalues_df) = c("Index", "Werte")
for (i in (0:4)){
  eigenvalues_df$Index[(i*M + 1): ((i+1)*M)] = rep((i+1), M)
  eigenvalues_df$Werte[(i*M + 1): ((i+1)*M)] = storage_ev[,(i+1)]
}
eigenvalues_df$Index <- as.factor(eigenvalues_df$Index)
# Plot
ggplot(eigenvalues_df, aes(x=Index, y=Werte)) + 
  geom_boxplot(color="blue",fill="blue", alpha=0.2, 
               notch=TRUE, notchwidth = 0.6,
               outlier.colour="red", outlier.fill="red", outlier.size=2) +
  annotate("segment", x = 0.5, xend = 4.5, y = 0.9, yend = 0.9, color = "darkgreen", linetype = "dashed", size = 0.8) +
  annotate("segment", x = 4.5, xend = 5.5, y = 1.4, yend = 1.4, color = "darkgreen", linetype = "dashed", size = 0.8) +
  labs(
    title = "Boxplot der Eigenwerte",
    x = "Index",
    y = "Eigenwerte"
  ) +
  theme_minimal()

# Histogram aller Eigenwerte
ggplot(eigenvalues_df, aes(x = Werte)) +
  geom_histogram(fill = "blue", color = "black", binwidth = 0.05, alpha = 0.2) + 
  geom_vline(aes(xintercept = mean(true_eigenvalues), color = "True Mean"), linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = mean(calc_ev), color = "Calculated Mean"), linetype = "solid", linewidth = 1) +
  scale_color_manual(values = c("True Mean" = "blue", "Calculated Mean" = "purple"), name = "Legende") +
  labs(
    title = "Histogram der Eigenwerte",
    x = "Werte",
    y = "Anzahl"
  ) + 
  theme_minimal()

# Histogram der Summen von Eigenwerten
sum_true_ev <- sum(true_eigenvalues)
sum_calc_ev <- rowSums(storage_ev)

ggplot() + 
  geom_histogram(aes(x = sum_calc_ev), binwidth = 0.05, color="black",fill="blue", alpha=0.2) +
  labs(x = "Summe der Eigenwerte", y = "Anzahl",
       title = "Histogram der Summen von Eigenwerten") +
  geom_vline(aes(xintercept = sum(true_eigenvalues), color = "True Sum"), linetype = "solid", linewidth = 1) +
  scale_color_manual(values = c("True Sum" = "blue"), name = "Legende") +
  theme_minimal()

