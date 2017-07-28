# Set file path
setwd("D:/F/Programming/Anomaly-Detection/Powergrid")

# Including function
if(!exists("formatMhsmm", mode="function")) source("../formatMHSMM.R")

# Read in data
train <- read.csv("train2008.csv", header = TRUE)
valid <- read.csv("train2009.csv", header = TRUE)
validAnom <- read.csv("valid20090.csv", header = TRUE)

VALIDATION_SIZE = 20000

# Extract the global_active_power column
trainGap <- formatMhsmm(data.frame(X = train$Global_active_power[1:20000]))
validGap <- formatMhsmm(data.frame(X = valid$Global_active_power[1:VALIDATION_SIZE]))
validAnomGap <- formatMhsmm(data.frame(X = validAnom$Noise[1:VALIDATION_SIZE]))
validAnomStatus<- formatMhsmm(data.frame(X = validAnom$Anomaly[1:VALIDATION_SIZE]))
# Load package
library(mvtnorm)
library(mhsmm)

# Build hmm with univariant data
# Set state number based on 10000 data in 

J <- 12
INTERVAL_SIZE = floor(VALIDATION_SIZE/J)

# Set init value pi
pi <- c(0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0834, 0.0834, 0.0834, 0.0834)

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

muPreset = c(sapply(data.frame(validAnom$Noise[1:INTERVAL_SIZE]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE+1:INTERVAL_SIZE*2]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*2+1:INTERVAL_SIZE*3]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*3+1:INTERVAL_SIZE*4]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*4+1:INTERVAL_SIZE*5]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*5+1:INTERVAL_SIZE*6]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*6+1:INTERVAL_SIZE*7]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*7+1:INTERVAL_SIZE*8]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*8+1:INTERVAL_SIZE*9]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*9+1:INTERVAL_SIZE*10]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*10+1:INTERVAL_SIZE*11]),mean, na.rm=TRUE),
             sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*11+1:INTERVAL_SIZE*12]),mean, na.rm=TRUE))

sigmaPreset = c(sapply(data.frame(validAnom$Noise[1:INTERVAL_SIZE]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE+1:INTERVAL_SIZE*2]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*2+1:INTERVAL_SIZE*3]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*3+1:INTERVAL_SIZE*4]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*4+1:INTERVAL_SIZE*5]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*5+1:INTERVAL_SIZE*6]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*6+1:INTERVAL_SIZE*7]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*7+1:INTERVAL_SIZE*8]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*8+1:INTERVAL_SIZE*9]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*9+1:INTERVAL_SIZE*10]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*10+1:INTERVAL_SIZE*11]),sd),
                sapply(data.frame(validAnom$Noise[INTERVAL_SIZE*11+1:INTERVAL_SIZE*12]),sd))

# Build emission matrix
B <- list(mu = muPreset,
          sigma = sigmaPreset)

# Build hmm model
model <- hmmspec(init = pi, trans = A, parms.emis = B, dens.emis = dnorm.hsmm)

# EM algorithom fits an HMM to the data
# @model: the init parameter to build hmm
hmm <- hmmfit(trainGap$x, model, mstep = mstep.norm, maxit = 150)

# Summary
summary(hmm)

# Test hmm model
yhat <- predict.hmm(hmm, validAnomGap$x, method = "viterbi")


# Define threshold value
threshold = 0.5
# build empty matrix
result_matrix <- matrix(NA, nrow = length(yhat$s), ncol = 4)

anomalies_correct = 0
anomalies_actual = 0
anomalies_found = 0
false_positives = 0
precision = 0
recall = 0

len <- length(yhat$s)-1

# point anomaly detection
for (i in 1:len) {
  mu = hmm$model$parms.emission$mu[yhat$s[i]]
  value <- abs(mu - validAnomGap$x[i])
  anomalous_value = validAnomGap$x[i]
  predicted_value = hmm$model$parms.emission$mu[yhat$s[i]]
  anomalous = validAnomStatus$x[i]
  if(anomalous == 1){
    anomalies_actual = anomalies_actual+1
  }
  if(value < threshold) {
    result_matrix[i,1] = 0
  } else{
    result_matrix[i,1] = 1
    anomalies_found = anomalies_found+1
    if(anomalous == 0){
      false_positives = false_positives+1
    }
    else{
      anomalies_correct = anomalies_correct+1
    }
  }
  result_matrix[i,2] <- value
  result_matrix[i,3] <- anomalous_value
  result_matrix[i,4] <- predicted_value
}

precision = anomalies_correct/anomalies_found
recall = anomalies_correct/anomalies_actual

# Write result to text file
write.table(result_matrix, file=sprintf("anomaly_%s.txt", threshold), 
            row.names = FALSE, col.names = FALSE)






