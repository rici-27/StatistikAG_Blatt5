library(tidyverse)
library(MASS)

## Aufgabe e)

Stock_Bond <- read_csv("Stock_Bond_2004_to_2006.csv", 
    col_types = cols(Date = col_skip(), DATE = col_date(format = "%m/%d/%Y"),
                     Three_month_treasury = col_skip(),
                     GM_Volume = col_skip(), F_Volume = col_skip(),
                     UTX_Volume = col_skip(), CAT_Volume = col_skip(),
                     MRK_Volume = col_skip(), PFE_Volume = col_skip(),
                     IBM_Volume = col_skip(), C_Volume = col_skip(), 
                     XOM_Volume = col_skip(), "S&P_Volume" = col_skip(),
                     MSFT_Volume = col_skip(),
                     "1 year Treasury Constant Maturity Rate" = col_skip(),
                     "3-Year Treasury Constant Maturity Rate" = col_skip(),
                     "10 year Treasury Constant Maturity Rate" = col_skip(),
                     "30 year Treasury Constant Maturity Rate"= col_skip(),
                     "Aaa Bond Yield" = col_skip(),                       
                     "Baa Bond Yield" = col_skip(),                        
                     "$/Euro" = col_skip(),                               
                     "Yen/$" = col_skip(),                                  
                     "Brazil Real/$" = col_skip()
                     ))
Data <- data.frame(Stock_Bond$DATE)
colnames(Data) <- c("Date")
Data$GM_returns <- c(NA, diff(log(Stock_Bond$GM_AC)))
Data$F_returns <- c(NA, diff(log(Stock_Bond$F_AC)))
Data$UTX_returns <- c(NA, diff(log(Stock_Bond$UTX_AC)))
Data$CAT_returns <- c(NA, diff(log(Stock_Bond$CAT_AC)))
Data$MRK_returns <- c(NA, diff(log(Stock_Bond$MRK_AC)))
Data$PFE_returns <- c(NA, diff(log(Stock_Bond$PFE_AC)))
Data$IBM_returns <- c(NA, diff(log(Stock_Bond$IBM_AC)))
Data$MSFT_returns <- c(NA, diff(log(Stock_Bond$MSFT_AC)))
Data$C_returns <- c(NA, diff(log(Stock_Bond$C_AC)))
Data$XOM_returns <- c(NA, diff(log(Stock_Bond$XOM_AC)))
Data$SP_returns <- c(NA, diff(log(Stock_Bond$SP_AC)))
# Zeile mit NA entfernen
Data <- Data[-1,]

# Data Frame mit Returns ohne Datum
Return <- Data[,-1]

# Wir betrachte AC 1-9
Return <- Return[,(1:9)]

# Gewichte berechnen
p <- ncol(Return)
days <- 336
n <- 336

# Funktion definieren, die Gewichte zurückgibt

get_weights <- function(X_transposed){
  # Gewichte Sample Cov berechnen
  sample_cov <- (1/n) * t(X_transposed) %*% X_transposed
  ones_vector <- c(rep(1,p))
  sample_cov_inverse <- solve(sample_cov)
  # sum() repräsentiert Skalarprodukt mit 1-Vektor im Nenner
  weights_sample_cov <- (sample_cov_inverse %*% ones_vector) / sum(sample_cov_inverse %*% ones_vector)
  
  # Gewichte Non Oracle berechnen
  # Parameter schätzen
  gamma <- (1/p) * sum(diag(sample_cov %*% diag(p)))
  delta2 <- (1/p) * norm((sample_cov - gamma * diag(p)), type = "F")^2
  beta2_help <- rep(0, n)
  for (i in 1:n){
    temp <- outer(t(X_transposed)[,i] , X_transposed[i,]) - sample_cov
    beta2_help[i] <- (1/p) * norm(temp, type = "F")^2 
  }
  beta2_pilot <- (1/n^2) * sum(beta2_help)
  beta2 <- min(beta2_pilot, delta2)
  alpha2 <- delta2 - beta2
  
  # Parameter berechnen
  para1 <- beta2/delta2
  para2 <- alpha2/delta2
  
  # Schätzer berechen
  non_oracle_est <- para1 * gamma * diag(p) + para2 * sample_cov
  non_oracle_inv <- solve(non_oracle_est)
  
  weights_non_oracle <- (non_oracle_inv %*% ones_vector) / sum(non_oracle_inv %*% ones_vector)
  
  weights_benchmark <- rep(1/p, p)
  
  return(list(weights_sample_cov = weights_sample_cov, 
              weights_non_oracle = weights_non_oracle,
              weights_benchmark = weights_benchmark))
}


# Speicher für Gewichte
# 1 für Sample Cov, 2 für Non Oracle, 3 für Benchmark
weights_storage <- array(NA, dim = c(days, 3, p))

for (i in (1:days)){
  X_transposed <- data.matrix(Return[i:(days+i-1),])
  weights_list <- get_weights(X_transposed)
  weights_storage[i,1,] <- weights_list$weights_sample_cov
  weights_storage[i,2,] <- weights_list$weights_non_oracle
  weights_storage[i,3,] <- weights_list$weights_benchmark
  print(i)
}


# Returns berechnen & vergleichen

return_sample_cov <- array(NA, dim = days+1)
return_non_oracle <- array(NA, dim = days+1)
return_benchmark <- array(NA, dim = days+1)
return_sample_cov[1] <- 0
return_non_oracle[1] <- 0
return_benchmark[1] <- 0

for (i in 1:days){
  return_sample_cov[i+1] <- return_sample_cov[i] +
                              sum(weights_storage[i,1, ] * Return[(days+i),] )
  return_non_oracle[i+1] <- return_non_oracle[i] +
                              sum(weights_storage[i,2, ] * Return[(days+i),] )
  return_benchmark[i+1] <- return_benchmark[i] +
                              sum(weights_storage[i,3, ] * Return[(days+i),] )
}

# Plot der Verläufe der Log Returns
ggplot() +
  geom_line(aes(x = Data$Date[(days):(2*days)], y = return_sample_cov, 
                color = "SampleCov"), size = 1) + 
  geom_line(aes(x = Data$Date[(days):(2*days)], y = return_non_oracle, 
                color = "NonOracle"), size = 1) + 
  geom_line(aes(x = Data$Date[(days):(2*days)], y = return_benchmark, 
                color = "Benchmark"), size = 1) + 
  labs(
    title = "Entwicklung der Returns unserer 3 Portfolien",
    x = "Datum",
    y = "Total Log Return"
  ) + 
  scale_color_manual(
    name = "Legende",
    values = c("SampleCov" = "blue",
               "NonOracle" = "purple",
               "Benchmark" = "black")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# Overall Return berechnen 

cat(
  "Overall Return für Sample Cov: ", return_sample_cov[days+1],"\n",
  "Overall Return für Non Oracle: ", return_non_oracle[days+1],"\n",
  "Overall Return für Benchmark: ", return_benchmark[days+1],"\n"
)

# Erwünschte Platzierung für AC 1-9 erreicht
