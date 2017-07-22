# Set file path
setwd("~/Desktop/Anomaly-Detection/Powergrid/")

# Read in data
train <- read.csv("train.csv", header = TRUE)
test <- read.csv("test_v1.csv")

# Extract the global_active_power column
Gap.col <- data.frame(X = train$Global_active_power[1:10000])
#Tgap.col <- data.frame(test$Global_active_power[1:20000])
# trainGap <- formatMhsmm(data.frame(X = train$Global_active_power[1:10000]))
# testGap <- formatMhsmm(data.frame(test$Global_active_power[764437:774437]))

# Replace 0 with NA
Gap.col[is.na(Gap.col)] <- 0

# Change data type
N <- as.numeric(Gap.col$X)
#result <- list(y = Gap.col$X)
class(N) <- "hsmm.data"

# Load package
library(mvtnorm)
library(mhsmm)

# Build hmm with univariant data

# Set state number based on 10000 data in graph
J <- 10

# Set init value pi
# #pi <- rep(1/J, J)
# pi <- c(0.111,0.050,0.081
#         ,0.090,0.075,0.035
#         ,0.250,0.170,0.138)
# pi <- c(runif(J, 0, 1))
# x <- replicate(J, diff(c(0, sort(runif(1)), 1)))
pi <- diff(c(0, sort(runif((J-1),0,1)), 1))

# Build transition matrix A, entry by row
# A <- matrix(c(0.8, 0.5, 0.1, 0.05, 0.2, 0.5, 0.15, 0.3, 0.4,
#               0.8, 0.5, 0.1, 0.05, 0.2, 0.5, 0.15, 0.4, 0.3,
#               0.8, 0.5, 0.1, 0.05, 0.4, 0.5, 0.15, 0.3, 0.2,
#               0.4, 0.5, 0.1, 0.05, 0.2, 0.5, 0.15, 0.3, 0.8,
#               0.8, 0.1, 0.5, 0.05, 0.2, 0.5, 0.15, 0.3, 0.4,
#               0.8, 0.5, 0.1, 0.05, 0.2, 0.15, 0.5, 0.3, 0.4,
#               0.15, 0.5, 0.1, 0.05, 0.2, 0.5, 0.8, 0.3, 0.4,
#               0.8, 0.05, 0.1, 0.5, 0.2, 0.5, 0.15, 0.3, 0.4,
#               0.8, 0.5, 0.1, 0.05, 0.5, 0.2, 0.15, 0.3, 0.4), nrow = J)

# A <- matrix(runif(J*J, 0,1),nrow = J)
A <- matrix(replicate(J,diff(c(0, sort(runif((J-1),0,1)), 1))),nrow = J)
base::t(A)   #transpose matrix
# Build emission matrix
B <- list(mu = c(runif(J, 0,1.5)), 
          sigma = c(runif(J, 0,0.5)))

# Build hmm model
model <- hmmspec(init = pi, trans = A, parms.emis = B, dens.emis = dnorm.hsmm, mstep = 10)
#model

# EM algorithom fits an HMM to the data
hmm <- hmmfit(N, model, mstep = mstep.norm, maxit = 500)

# Summary
summary(hmm)

plot(hmm$loglik, type = 'l', ylab = "log-likelihood", xlab = "Iteration")
yhat <- predict(hmm, N)
plot(yhat)



# Test hmm model
# yhat1 <- predict(hmm, trainGap$x)
# yhat2 <- predict(hmm, testGap$x)
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



