# Set file path
setwd("~/Desktop/Anomaly-Detection/Powergrid/")

# Including function
if(!exists("formatMhsmm", mode="function")) source("../formatMHSMM.R")

# Read in data
train <- read.csv("train.csv", header = TRUE)
test <- read.csv("test_v1.csv", header = TRUE)

# Extract the global_active_power column
trainGap <- formatMhsmm(data.frame(X = train$Global_active_power[1:300000]))
testGap <- formatMhsmm(data.frame(test$Global_active_power[764437:774437]))

# Load package
library(mvtnorm)
library(mhsmm)

# Build hmm with univariant data

# Set state number based on 10000 data in graph
J <- 15

# Set init value pi
#pi <- rep(1/J, J)
#pi <- c(0.1, 0.75, 0.15)
pi <- c(0.123, 0.107, 0.088, 0.036, 0.099,
        0.001, 0.025, 0.095, 0.006, 0.018,
        0.085, 0.092, 0.024, 0.120, 0.081)


# Build transition matrix A, entry by row
# A <- matrix(c(0.9, 0.5, 0.1, 
#               0.08, 0.2, 0.5, 
#               0.02, 0.3, 0.4), nrow = J)

A <- matrix(c(0.104,	0.13,	 0.014,	0.115,	0.041,	0.103,	0.029,	0.052,	0.106,	0.12,	  0.046,	0.14,	  0.105,	0.053,	0.021,
              0.119,	0.096, 0.005,	0.012,	0.104,	0.053,	0.137,	0.099,	0.018,	0.105,	0.032,	0.136,	0.08,	  0.173,	0.097,
              0.028,	0.018, 0.041,	0.092,	0.103,	0.062,	0.036,	0.101,	0.072,	0.12,	  0.002,	0.022,	0.036,  0.065,	0.024,
              0.047,	0.015, 0.039,	0.014,	0.058,	0.099,	0.039,	0.032,	0.12,	  0.005,	0.113,	0.115,	0.059,	0.039,	0.081,
              0.12,	  0.137, 0.003,	0.019,	0.048,	0.043,	0.127,	0.046,	0.042,	0.042,	0.115,	0.048,	0.086,	0.053,	0.032,
              0.058,	0.019, 0.002,	0.063,	0.089,	0.055,	0.055,	0.048,  0.04,	  0.069,	0.001,	0.053,	0.05,	  0.072,	0.056,
              0.032,	0.03,	 0.18,	0.059,	0.047,	0.011,	0.014,	0.101,	0.071,	0.085,	0.05,	  0.115,	0.08,	  0.065,	0.118,
              0.011,	0.03,	 0.061,	0.036,	0.103,	0.056,	0.049,	0.049,	0.083,	0.032,	0.129,	0.029,	0.051,	0.012,	0.042,
              0.026,	0.08,	 0.205,	0.035,	0.032,	0.029,	0.084,	0.01,	  0.062,	0.051,	0.055,	0.055,	0.051,	0.051,	0.094,
              0.085,	0.146, 0.036,	0.064,	0.133,	0.012,  0.067,	0.105,	0.072,	0.019,	0.067,	0.028,	0.015,	0.176,	0.036,
              0.067,	0.081, 0.126,	0.1,	  0.029,	0.134,	0.09,	  0.032,	0.049,	0.11,	  0.015,	0.103,	0.119,	0.028,	0.065,
              0.058,	0.122, 0.028,	0.125,	0.044,	0.043,	0.017,	0.063,	0.005,	0.052,	0.115,	0.004,	0.063,	0.09,	  0.03,
              0.078,	0.074, 0.209,	0.102,	0.036,	0.143,	0.039,	0.051,	0.123,	0.114,	0.049,	0.039,	0.045,	0.021,	0.097,
              0.121,	0.002, 0.043,	0.082,	0.06,	  0.084,	0.148,	0.096,	0.012,	0.048,	0.082,	0.033,	0.095,	0.058,	0.116,
              0.046,	0.02,	 0.008,	0.082,	0.073,	0.073,	0.069,	0.115,	0.125,	0.028,	0.129,	0.08,	  0.065,	0.044,	0.091), nrow = J)

# Build emission matrix
#B <- list(mu = c(1.26, 1.86, 2.6), sigma = c(0, 4, 8))

B <- list(mu = c(0.16, 0.23, 0.283, 0.329, 0.391, 
                 0.457, 0.51, 0.539, 0.57, 0.622,
                 0.734, 0.927, 1.184, 1.343, 1.456), 
          sigma = c(0.043, 0.011, 0.017, 0.013, 0.017,
                    0.021, 0.01, 0.008, 0.01, 0.022, 
                    0.043, 0.069,0.07, 0.031, 0.038))

# Build hmm model
model <- hmmspec(init = pi, trans = A, parms.emis = B, dens.emis = dnorm.hsmm, mstep = 10)
#model

# EM algorithom fits an HMM to the data
hmm_train <- hmmfit(trainGap$x, model, mstep = mstep.norm, maxit = 500)
#hmm_test <- hmmfit(testGap$x, model, mstep = mstep.norm, maxit = 500)

# Summary
summary(hmm_train)
plot(hmm_train$loglik, col = "red", type = 'b', 
     ylab = "log-likelihood", xlab = "Iteration")
#points(hmm_test$loglik, col = "blue", type = 'b')

# Test hmm model
yhat1 <- predict(hmm_train, trainGap$x)
yhat2 <- predict(hmm_train, testGap$x)

# Plot
plot(yhat1, type = 'l', ylab ='Global_active_power_train')
plot(yhat2, type = 'l', ylab ='Global_active_power_test')

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



