# Set file path
path<-getwd()
setwd(path)

# Including function
if(!exists("formatMhsmm", mode="function")) source("../formatMHSMM.R")

# Read in data
train <- read.csv("train.csv", header = TRUE)
vali<-read.csv("vali2.csv", header = TRUE)

# Extract the global_active_power column
trainGap <- formatMhsmm(data.frame(X = train$Global_active_power[1:524556]))
valiRealGap <- formatMhsmm(data.frame(X = train$Global_active_power[524557:569196]))
valiAnomGap<- formatMhsmm(data.frame(X = vali$Noise[1:44640]))
valiAnomStatus<- formatMhsmm(data.frame(X = vali$Anomaly[1:44640]))

# Load package
library(mvtnorm)
library(mhsmm)

# Build hmm with univariant data
# Set state number based on 10000 data in graph
J <- 6

# Set init value pi
pi <- c(0.5, 0.2, 0.1, 0.1, 0.08, 0.02)

A <- matrix(c(0.6,	0.25,	0.1,	0.01,	0.01,	0.32,
              0.2,	0.5,	0.16,	0.02,	0.01,	0.02,
              0.15,	0.2,	0.5,	0.13,	0.01,	0.03,
              0.02,	0.03,	0.21,	0.5,	0.1,	0.05,
              0.02,	0.01,	0.01,	0.21,	0.55,	0.16,
              0.01,	0.01,	0.02,	0.13,	0.32,	0.42), nrow = J)

# Build emission matrix
B <- list(mu = c(1, 3, 5, 7, 9, 11), 
          sigma = c(1,1,1,1,1,1))

# Build hmm model
model <- hmmspec(init = pi, trans = A, parms.emis = B, dens.emis = dnorm.hsmm)

# EM algorithom fits an HMM to the data
# @model: the init parameter to build hmm
hmm <- hmmfit(trainGap$x, model, mstep = mstep.norm, maxit = 500)

# Summary
summary(hmm)

# Test hmm model
#yhat1 <- predict.hmm(hmm, trainGap$x, method = "viterbi")
yhat <- predict.hmm(hmm, valiRealGap$x, method = "viterbi")


# Define threshold value
threshold = 0.5
# build empty matrix
result_matrix <- matrix(NA, nrow = length(valiRealGap$x), ncol = 4)
anomalies_correct = 0
anomalies_actual = 0
anomalies_found = 0
false_positives = 0
precision = 0
recall = 0

# point anomaly detection
for (i in 1:length(valiRealGap$x)) {
  anomalous_value = valiAnomGap$x[i]
  predicted_value = hmm$model$parms.emission$mu[yhat$s[i]]
  anomalous = valiAnomStatus$x[i]
  value = round(abs(anomalous_value - predicted_value), 3)
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

