library(readr)

Stock_Bond <- read_csv("Stock_Bond_2004_to_2006.csv", 
    col_types = cols(Date = col_skip(), DATE = col_date(format = "%m/%d/%Y"),
                     Date = col_skip(), Three_month_treasury = col_skip(), GM_Volume = col_skip(),
                     F_Volume = col_skip(), UTX_Volume = col_skip(), CAT_Volume = col_skip(),
                     MRK_Volume = col_skip(), PFE_Volume = col_skip(), IBM_Volume = col_skip(),
                     C_Volume = col_skip(), XOM_Volume = col_skip(), "S&P_Volume" = col_skip(),
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
Data <- Data[-1,]
Data <- Data[,-1]

# Gewichte berechnen

p <- 11
days <- 336
n <- 336

# Funktionen definieren

get_sample_weights <- function(X){
  sample_cov <- (1/n) * t(X) %*% X
  ones_vector <- c(rep(1,p))
  sample_cov_inverse <- solve(sample_cov)
  weights <- (sample_cov_inverse %*% ones_vector) / sum(sample_cov_inverse %*% ones_vector)
  return(weights)
}

get_non_oracle_weights <- function(X){
  sample_cov <- (1/n) * t(X) %*% X

  # Parameter schätzen
  gamma <- (1/p) * sum(diag(sample_cov %*% diag(p)))
  delta2 <- (1/5) * norm((sample_cov - gamma * diag(p)), type = "F")^2
  beta2_help <- rep(0, n)
  for (i in 1:n){
    temp <- outer(t(X)[,i] , X[i,]) - sample_cov
    beta2_help[i] <- (1/p) * norm(temp, type = "F")^2 
  }
  beta2_pilot <- (1/n^2) * sum(beta2_help)
  beta2 <- min(beta2_pilot, delta2)
  alpha2 <- delta2 - beta2
  
  # Gewichte berechnen
  rho1 <- (beta2/delta2) * gamma
  rho2 <- alpha2/delta2
  
  # Schätzer berechen
  non_oracle_est <- rho1 * diag(p) + rho2 * sample_cov
  non_oracle_inv <- solve(non_oracle_est)
  
  #Gewichte
  ones_vector <- c(rep(1,p))
  weights <- (non_oracle_inv %*% ones_vector) / sum(non_oracle_inv %*% ones_vector)
  return(weights)
}

# Speicher für Gewichte
# 1 für Sample Cov, 2 für Non Oracle, 3 für Benchmark
weights_storage <- array(NA, dim = c(days, 3, p))
weights_storage[,3,] <- 1/p

for (i in (1:days)){
  X <- data.matrix(Data[i:(days+i-1),])
  weights_storage[i,1,] <- get_sample_weights(X)
  weights_storage[i,2,] <- get_non_oracle_weights(X)
}

View(weights_storage[1,,])

# Overall return berechnen & vergleichen

overall_return_sample_cov <- 0
overall_return_non_oracle <- 0
overall_return_benchmark <- 0

for (i in 1:days){
  overall_return_sample_cov <- overall_return_sample_cov + sum(weights_storage[i,1, ] * Data[(days+i),] )
  overall_return_non_oracle <- overall_return_shrinkage + sum(weights_storage[i,2, ] * Data[(days+i),] )
  overall_return_benchmark <- overall_return_benchmark + sum(weights_storage[i,3, ] * Data[(days+i),] )
}

print(paste("Overall return for sample cov: " ,  as.character(overall_return_sample_cov) , " ."))
print(paste("Overall return for shrinkage: " ,  as.character(overall_return_shrinkage) , " ."))
print(paste("Overall return for benchmark: " ,  as.character(overall_return_benchmark) , " ."))




