library(tidyverse) # für ggplot2 zum plotten
library(MASS) # Simulation von multivariat normalverteilten Daten

# Vorbereitung

n <- 100
p <- 5
# Kovarianz Matrix und Erwartungswert erstellen
Sigma <- diag(0.9, nrow = p, ncol = p) + matrix(0.1, nrow = p, ncol = p)
mu <- c(rep(0,p))
# Multivariate Daten erzeugen
X_transposed <- mvrnorm(n, mu, Sigma)
# X_transposed in R^{n \times p}


## Aufgabe a)

M <- 1000
# symmetric -> Matrix ist symmetrisch
# only.values -> Eigenvektoren werden nicht berechnet
# decreasing = FALSE -> Eigenwerte werden aufsteigend sortiert
true_eigenvalues <- sort(eigen(Sigma, symmetric=TRUE,
                               only.values=TRUE)$values, decreasing = FALSE)

# Speicher für Eigenwerte erstellen
storage_ev <- matrix(NA, nrow = M, ncol = p)

# Monte Carlo Simulation
for (i in (1:M)){
  X_transposed <- mvrnorm(n, mu, Sigma)
  S <- (1/n) * t(X_transposed) %*% X_transposed
  storage_ev[i,] <- sort(eigen(S, symmetric=TRUE,
                               only.values=TRUE)$values, decreasing = FALSE)
}

# Boxplot der Eigenwerte
# Zuerst Data Frame erstellen
eigenvalues_df <- data.frame(matrix(NA, nrow=p*M, ncol=2))
colnames(eigenvalues_df) = c("Index", "Werte")
for (i in (0:4)){
  eigenvalues_df$Index[(i*M + 1): ((i+1)*M)] = rep((i+1), M)
  eigenvalues_df$Werte[(i*M + 1): ((i+1)*M)] = storage_ev[,(i+1)]
}
# Index Spalte von numerischen in kategorische Werte umwandeln
eigenvalues_df$Index <- as.factor(eigenvalues_df$Index)
View(eigenvalues_df)

# Box Plot der Eigenwerte
ggplot(eigenvalues_df, aes(x = Index, y = Werte)) + 
  geom_boxplot(color = "blue", fill = "blue", alpha = 0.2, notch = FALSE,
               outlier.colour = "red", outlier.fill = "red", outlier.size = 2) +
  geom_segment(aes(x = 0.5, xend = 4.5, y = 0.9, yend = 0.9, color = "Eigenwert 0.9"), 
               linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE) +
  geom_segment(aes(x = 4.5, xend = 5.5, y = 1.4, yend = 1.4, color = "Eigenwert 1.4"), 
               linetype = "dashed", linewidth = 0.8, inherit.aes = FALSE) +
  # inherit.as = FALSE -> Wir überschreiben das globale aes
  scale_color_manual(
    values = c("Eigenwert 0.9" = "darkgreen", "Eigenwert 1.4" = "darkred"),
    name = "Legende"
  ) +
  labs(
    title = "Boxplot der Eigenwerte",
    x = "Index",
    y = "Eigenwerte"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  )

# Histogram der Summen von Eigenwerten
sum_true_ev <- sum(true_eigenvalues)
sum_calc_ev <- rowSums(storage_ev)

ggplot() + 
  geom_histogram(aes(x = sum_calc_ev), binwidth = 0.05,
                 color="black", fill="blue", alpha=0.2) +
  labs(x = "Summe der Eigenwerte", y = "Anzahl",
       title = "Histogramm der Summen der Eigenwerte") +
  geom_vline(aes(xintercept = sum_true_ev, color = "True Sum"),
             linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = mean(sum_calc_ev), color = "Mean Sum"),
             linetype = "dotted", linewidth = 1) +
  scale_color_manual(values = c("True Sum" = "blue", "Mean Sum" = "purple"), 
                     labels = c(
                       paste0("True Sum: ", sum_true_ev),
                       paste0("Mean Sum: ", round(mean(sum_calc_ev), 4))),
                     name = "Legende") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )
  
# Kleines Add-On: 
# Histogram aller Eigenwerte
ggplot(eigenvalues_df, aes(x = Werte)) +
  geom_histogram(fill = "blue", color = "black", binwidth = 0.05, alpha = 0.2) + 
  geom_vline(aes(xintercept = mean(true_eigenvalues), color = "True Mean"),
             linetype = "solid", linewidth = 1) +
  geom_vline(aes(xintercept = mean(storage_ev), color = "Calculated Mean"),
             linetype = "dotted", linewidth = 1) +
  scale_color_manual(values = c("True Mean" = "blue", "Calculated Mean" = "purple"),
                     labels = c(
                       paste0("True Mean: ", mean(true_eigenvalues)),
                       paste0("Mean: ", round(mean(storage_ev),4))),
                     name = "Legende") + 
  labs(
    title = "Histogramm der Eigenwerte",
    x = "Werte",
    y = "Anzahl"
  ) + 
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )
