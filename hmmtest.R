# Set file path
path<-getwd()
setwd(path)

# Read in data
train <- read.csv("train.csv", header = TRUE)
#test <- read.csv("test_v1.csv")

# Extract the global_active_power column
Gap.col <- data.frame(X = train$Global_active_power[134318:661358])
#Tgap.col <- data.frame(test$Global_active_power[1:20000])

NoNa <- na.omit(Gap.col$X)

# Change data type
N <- as.numeric(NoNa)
#result <- list(y = Gap.col$X)
class(N) <- "hsmm.data"

# Load package
library(mvtnorm)
library(mhsmm)

# Build hmm with univariant data

# Set state number based on 10000 data in graph
J <- 12

# Set init value pi+
#pi <- rep(1/J, J)
pi <- c(0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0834, 0.0834, 0.0834, 0.0834)

# Build transition matrix A, entry by row
A <- matrix(c(0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
              0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834,
              0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834,
              0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834,
              0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834, 0.0834), nrow = J)

# Build emission matrix
B <- list(mu = c(0.4850526, 1.4811457, 2.9129929,
                 0.3850296, 1.1706618, 2.3604874,
                 0.4556026, 1.5080793, 3.0992803,
                 0.4531069, 1.5419564, 3.2906781), 
          sigma = c(0.03678481, 0.15216245, 1.03251746,
                    0.03502272, 0.18842031, 0.76244307,
                    0.02878864, 0.15348746, 1.12532468,
                    0.02581237, 0.22541860, 1.26231614
                    ))

# Build hmm model
model <- hmmspec(init = pi, trans = A, parms.emis = B, dens.emis = dnorm.hsmm, mstep = 10)
model

# EM algorithom fits an HMM to the data
hmm <- hmmfit(N, model, mstep = mstep.norm, maxit = 200)

# Summary
summary(hmm)
print(hmm$loglik)

#plot(hmm$loglik, type = 'l', ylab = "log-likelihood", xlab = "Iteration")
yhat <- predict(hmm, N)
plot(yhat)


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



