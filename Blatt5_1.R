library(tidyverse)
library(MASS)

# Vorbereitung

n <- 100
dim <- 5
# Kovariance Matrix und Mittelwert erstellen
Sigma <- matrix(c( 1 , rep (.1 ,4) , .1 , 1 , rep ( .1 , 3 ) , .1 , .1 , 1 , .1 , .1 , rep(.1 ,3) ,1 ,.1 , rep (.1 ,4) ,1) ,5 ,5)
mu <- c(rep(0,dim))
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
storage_ev <- matrix(0, nrow = M, ncol = dim)
# Monte Carlo Simulation
for (i in (1:M)){
  X <- mvrnorm(n, mu, Sigma)
  S <- (1/n) * t(X) %*% X
  storage_ev[i,] <- sort(eigen(S)$values)
}
# Mittelwerte der Eigenwerte berechnen
calc_ev <- (1/M) * colSums(storage_ev)

# Einfacher Plot von berechneten und wahren Eigenwerten

indizes <- c(1,2,3,4,5)
ggplot() +
  geom_point(aes(x=indizes, y = true_eigenvalues), size = 2, color = "blue") +
  geom_point(aes(x=indizes, y = calc_ev), size = 2, color = "darkgreen") +
  labs(x = "Indizes", y = "Eigenwerte", title = "Wahre & Berechnete Eigenwerte") +
  theme_minimal()

# Boxplot
# Data Frame erstellen
### hier vllt noch wahre werte eintragen?
eigenvalues_df <- data.frame(matrix(0, nrow=dim*M, ncol=2))
colnames(eigenvalues_df) = c("Index", "Werte")
for (i in (0:4)){
  eigenvalues_df$Index[(i*M + 1): ((i+1)*M)] = rep((i+1), M)
  eigenvalues_df$Werte[(i*M + 1): ((i+1)*M)] = storage_ev[,(i+1)]
}
eigenvalues_df$Index <- as.factor(eigenvalues_df$Index)

ggplot(eigenvalues_df, aes(x=Index, y=Werte)) + 
  geom_boxplot(color="blue",fill="blue", alpha=0.2, 
               notch=TRUE, notchwidth = 0.6,
               outlier.colour="red", outlier.fill="red", outlier.size=2) +
  labs(
    title = "Boxplot der Eigenwerte",
    x = "Index",
    y = "Eigenwerte"
  ) +
  theme_minimal()

# Histogram der Eigenwerte 
# weg? legende hinzufügen!
ggplot(eigenvalues_df, aes(x=Werte)) +
  geom_histogram(fill="lightpink", binwidth = 0.05) +
  geom_vline(xintercept = mean(true_eigenvalues), color = "blue", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = mean(c(calc_ev)), color = "purple", linetype = "dotted", linewidth = 1) +
  labs(
    title = "Histogram der Eigenwerte",
    x = "Werte",
    y = "Anzahl"
  ) +
  theme_minimal()
# legen und überschrift fehlt usw.

sum_true_ev <- sum(true_eigenvalues)
sum_calc_ev <- rowSums(storage_ev)

ggplot() + 
  geom_histogram(aes(x = sum_calc_ev), binwidth = 0.05, color="blue",fill="blue", alpha=0.4) +
  labs(x = "Summe der Eigenwerte", y = "Anzahl",
       title = "Histogram der Summen von Eigenwerten") +
  theme_minimal()

# vllt auch sum in histogram packen
