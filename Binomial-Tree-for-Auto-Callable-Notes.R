
S <- 108.04
r <- 0.0223
r1 <- 0.0285
div <- 0.46/S
sigma <- 0.24542
t <- 3
threshold <- c(1.05*S, 1.1*S, 1.15*S, 0.8*S)
times <- 10
error = vector()
Value = vector()
determination1 <- c(94,
                    185,
                    273,
                    367,
                    458,
                    549,
                    640,
                    731,
                    823,
                    915,
                    1004)
pay1 <- c(97,
          188,
          278,
          370,
          461,
          552,
          643,
          734,
          826,
          920,
          1007)
divdate1 <- c(12,
              103,
              201,
              285,
              376,
              474,
              565,
              656,
              747,
              838,
              929,
              1020)

#Steps from 1*1096 to 10*1096
for (times in 1:10) {
  
  determination <- determination1 * times
  pay <- pay1 * times
  divdate <- divdate1 * times
  
  VStock <- matrix(0, nrow=1096*times+1, ncol=1096*times+1)
  VOption <- matrix(0, nrow=1096*times+1, ncol=1096*times+1)
  
  n <- 1096*times
  #parameters
  delta <- t/n
  u <- exp(sigma*sqrt(delta))
  d <- exp(-sigma*sqrt(delta))
  qu <- (exp(r1*delta)-d)/(u-d)
  qd <- 1-qu
  #stock tree 
  k <- 0
  for (j in 1:(n+1)) {
    if (is.element(j, divdate)) {k<-k+1}
    for (i in 1:j) {
      VStock[j,i]=S*u^(i-1)*d^(j-i)*(1-div)^k
    }
  }
  #option tree
  #maturity
  j <- n+1
  for (i in 1:j) {
    VOption[j,i]=ifelse(VStock[j,i]>=threshold[4], 10*(1+0.0895/4)*exp(-r*3/365), 10*VStock[j,i]/S*exp(-r*3/365))
  }
  #before maturity
  l <- 0
  for (j in n:1) {
    if (is.element(j, determination)) {
      l <- which(determination==j)
      B <- threshold[(l-1)%/%4+1]
      for (i in 1:j) {
        if (VStock[j,i] < B) {
          if (VStock[j,i] >= threshold[4]) {
            VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])+10*0.0895/4*exp(-r*(pay[l]-determination[l])/365)
          } else {
            VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])
          }
        } else {
          VOption[j,i]=10*(1+0.0895/4)*exp(-r*(pay[l]-determination[l])/365)
        }
      }
    } else {
      for (i in 1:j) {
        VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])
      }
    }
  }
  
  Value[times] = VOption[1,1]
  
  #error
  error[times] <- VOption[1,1] - 9.646

}

print(Value)

rm(list=ls())



#Sigma Sensitivity Analysis

S <- 108.04
r <- 0.0223
r1 <- 0.0285
div <- 0.46/S
sigma <- 0.2
t <- 3
threshold <- c(1.05*S, 1.1*S, 1.15*S, 0.8*S)
error = vector()
Value = vector()
determination <- c(94,
                   185,
                   273,
                   367,
                   458,
                   549,
                   640,
                   731,
                   823,
                   915,
                   1004)
pay <- c(97,
         188,
         278,
         370,
         461,
         552,
         643,
         734,
         826,
         920,
         1007)
divdate <- c(12,
             103,
             201,
             285,
             376,
             474,
             565,
             656,
             747,
             838,
             929,
             1020)
times <- 1

count <- 0

sigmasensitivity = vector()

for (sigma in seq(0.2, 0.3, 0.001)) {
  count <- count + 1
  
  VStock <- matrix(0, nrow=1096*times+1, ncol=1096*times+1)
  VOption <- matrix(0, nrow=1096*times+1, ncol=1096*times+1)
  
  n <- 1096*times
  
  #parameters
  delta <- t/n
  u <- exp(sigma*sqrt(delta))
  d <- exp(-sigma*sqrt(delta))
  qu <- (exp(r1*delta)-d)/(u-d)
  qd <- 1-qu
  
  #stock tree 
  k <- 0
  for (j in 1:(n+1)) {
    if (is.element(j, divdate)) {k<-k+1}
    for (i in 1:j) {
      VStock[j,i]=S*u^(i-1)*d^(j-i)*(1-div)^k
    }
  }
  
  #option tree
  #maturity
  j <- n+1
  for (i in 1:j) {
    VOption[j,i]=ifelse(VStock[j,i]>=threshold[4], 10*(1+0.0895/4)*exp(-r*3/365), 10*VStock[j,i]/S*exp(-r*3/365))
  }
  
  #before maturity
  l <- 0
  for (j in n:1) {
    if (is.element(j, determination)) {
      l <- which(determination==j)
      B <- threshold[(l-1)%/%4+1]
      for (i in 1:j) {
        if (VStock[j,i] < B) {
          if (VStock[j,i] >= threshold[4]) {
            VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])+10*0.0895/4*exp(-r*(pay[l]-determination[l])/365)
          } else {
            VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])
          }
        } else {
          VOption[j,i]=10*(1+0.0895/4)*exp(-r*(pay[l]-determination[l])/365)
        }
      }
    } else {
      for (i in 1:j) {
        VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])
      }
    }
  }
  
  sigmasensitivity[count] <- VOption[1,1]
  
}  

plot(x=seq(0.2, 0.3, 0.001), y=sigmasensitivity, type = 'l', xlab='sigma', ylab='Option Value')
abline(h=9.646, col="red", lty = 2)
legend("topright", c("Option Value", "Estimated Value"), col = c("black", "red"), lty = c(1,2))

rm(list=ls())



#Dividend Sensitivity Analysis

S <- 108.04
r <- 0.0223
r1 <- 0.0285
div <- 0.2/S
sigma <- 0.24542
t <- 3
threshold <- c(1.05*S, 1.1*S, 1.15*S, 0.8*S)
error = vector()
Value = vector()
determination <- c(94,
                   185,
                   273,
                   367,
                   458,
                   549,
                   640,
                   731,
                   823,
                   915,
                   1004)
pay <- c(97,
         188,
         278,
         370,
         461,
         552,
         643,
         734,
         826,
         920,
         1007)
divdate <- c(12,
             103,
             201,
             285,
             376,
             474,
             565,
             656,
             747,
             838,
             929,
             1020)

times <- 1
count <- 0
dividend <- 0.2


divsensitivity = vector()

for (dividend in seq(0.2, 0.6, 0.01)) {
  count <- count + 1
  
  VStock <- matrix(0, nrow=1096*times+1, ncol=1096*times+1)
  VOption <- matrix(0, nrow=1096*times+1, ncol=1096*times+1)
  
  n <- 1096*times
  
  #parameters
  delta <- t/n
  u <- exp(sigma*sqrt(delta))
  d <- exp(-sigma*sqrt(delta))
  qu <- (exp(r1*delta)-d)/(u-d)
  qd <- 1-qu
  
  #stock tree 
  k <- 0
  for (j in 1:(n+1)) {
    if (is.element(j, divdate)) {k<-k+1}
    for (i in 1:j) {
      VStock[j,i]=S*u^(i-1)*d^(j-i)*(1-div)^k
    }
  }
  
  #option tree
  #maturity
  j <- n+1
  for (i in 1:j) {
    VOption[j,i]=ifelse(VStock[j,i]>=threshold[4], 10*(1+0.0895/4)*exp(-r*3/365), 10*VStock[j,i]/S*exp(-r*3/365))
  }
  
  #before maturity
  l <- 0
  for (j in n:1) {
    if (is.element(j, determination)) {
      l <- which(determination==j)
      B <- threshold[(l-1)%/%4+1]
      for (i in 1:j) {
        if (VStock[j,i] < B) {
          if (VStock[j,i] >= threshold[4]) {
            VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])+10*0.0895/4*exp(-r*(pay[l]-determination[l])/365)
          } else {
            VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])
          }
        } else {
          VOption[j,i]=10*(1+0.0895/4)*exp(-r*(pay[l]-determination[l])/365)
        }
      }
    } else {
      for (i in 1:j) {
        VOption[j,i]=exp(-r*delta)*(qu*VOption[j+1,i+1]+qd*VOption[j+1,i])
      }
    }
  }
  
  divsensitivity[count] <- VOption[1,1]
  dividend <- dividend + 0.01
  div <- dividend/S
  
}  

plot(x=seq(0.2, 0.6, 0.01), y=divsensitivity, type = 'l', xlab='Dividend', ylab='Option Value')
abline(h=9.646, col="red", lty = 2)
legend("topright", c("Option Value", "Estimated Value"), col = c("black", "red"), lty = c(1,2))


rm(list=ls())
