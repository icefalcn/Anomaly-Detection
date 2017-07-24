# Set file path
setwd("~/Desktop/Anomaly-Detection/Powergrid/")

# Including function
if(!exists("formatMhsmm", mode="function")) source("../formatMHSMM.R")

# Read in data
train <- read.csv("train.csv", header = TRUE)
test <- read.csv("test_v1.csv", header = TRUE)

# Extract the global_active_power column
trainGap <- formatMhsmm(data.frame(X = train$Global_active_power))
testGap <- formatMhsmm(data.frame(test$Global_active_power))

# Load package
library(mvtnorm)
library(mhsmm)

# Build hmm with univariant data
# Set state number based on 10000 data in graph
J <- 6

# Set init value pi
pi <- c(1, 0, 0, 0, 0, 0)

A <- matrix(c(0.6,	0.25,	0.1,	0.01,	0.01,	0.32,
              0.2,	0.5,	0.16,	0.02,	0.01,	0.02,
              0.15,	0.2,	0.5,	0.13,	0.01,	0.03,
              0.02,	0.03,	0.21,	0.5,	0.1,	0.05,
              0.02,	0.01,	0.01,	0.21,	0.55,	0.16,
              0.01,	0.01,	0.02,	0.13,	0.32,	0.42), nrow = J)

# Build emission matrix
B <- list(mu = c(1.204, 1.02, 0.754, 1.334, 1.731, 1.514), 
          sigma = c(0.989, 0.84, 0.785, 1.14, 1.302, 1.279))

# Build hmm model
modelspec <- hmmspec(init = pi, trans = A, parms.emis = B, dens.emis = dnorm.hsmm, mstep = 10)

# EM algorithom fits an HMM to the data
hmm_model <- hmmfit(trainGap$x, modelspec, mstep = mstep.norm, maxit = 500)

# Summary
summary(hmm_model)

# Test hmm model
yhat1 <- predict.hmm(hmm_model, trainGap$x)
yhat2 <- predict.hmm(hmm_model, testGap$x)

# plot(yhat1$x, col = 'red', type = 'b', ylab ='Global_active_power_train')
# points(yhat2$x, col = 'cyan', type = 'b')
# plot(yhat2, type = 'l', ylab ='Global_active_power_test')

plot(yhat1$x[1:20000], col = 'red', type = 'l')
points(yhat2$x[1:10000], col = 'blue', type = 'p')
legend('topright', legend=c("data by HMM", "observation"), 
       col=c('red', 'blue'), lty = c(1, NA), pch = c(NA, 1), cex = 0.6)

# Define threshold value
threshold = 0.5
# build empty matrix
result_matrix <- matrix(NA, nrow = length(testGap$x), ncol = 2)

# point anomaly detection
for (i in 1:length(testGap$x)) {
  value = round(abs(yhat1$x[i] - yhat2$x[i]), 3)
      if(value < threshold) {
        result_matrix[i,1] = 0
      } else{
        result_matrix[i,1] = 1
      }
      result_matrix[i,2] <- value
}

# Write result to text file
write.table(result_matrix, file=sprintf("anomaly_%s.txt", threshold), 
            row.names = FALSE, col.names = FALSE)

# J <- 3
# init <- c(0,0,1)
# P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
# B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
# d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
# model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
# train <- simulate(x, model,r=rnorm.hsmm,nsim=100,seed=123456)
# plot(train,xlim=c(0,400))
# 
# predict()



