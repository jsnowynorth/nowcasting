####################################################################
### Joshua North
### 10/15/2018
### MSE forecast horizon
####################################################################

####################################################################
### Load Libraries
####################################################################
library(fields)
library(ncdf4)
library(RNetCDF)
library(xtable)
library(forecast)
library(plot3D)
library(ggplot2)
library(MASS)
library(mvtnorm)

####################################################################
### Souce Functions
####################################################################
source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Load_netCDF_Data.R', chdir = TRUE)
source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Nowcasting_Function.R', chdir = TRUE)
source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Finite_Difference_Function.R', chdir = TRUE)
source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Move_Functions.R', chdir = TRUE)
source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Error_Function.R', chdir = TRUE)
source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Prediction_Functions.R', chdir = TRUE)


####################################################################
### Compile Data
####################################################################

bigD <- load.netCDF.data(filepath = '/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/lightning_data')

run <- read.csv('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Movement_Data.csv') # mac
run <- run[,2:5] # drop first column

move.perc <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/predict_data.csv', header = TRUE)
col.prec <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/column_predict.csv', header = TRUE)


Q <- (bigD[,,60])/100
t1_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t1 <- finiteD(deltaT = .01, deltaX = .0025, uInit = t1_a$obs1, nframes = 20)



###### MSE Estimates

for(j in 3:62){
  compNQ <- matrix(c(0), nc = 5, nr = 10)
  colnames(compNQ) <- c("Timestep", "FDAWLM", "Persisted", "Difference", "Change")
  
  begin = j
  q <- (bigD[,,begin] - bigD[,,(begin-1)])/100
  #Q <- ifelse(Q < 0.15, 0, Q)
  sampsNQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                     timestep = begin, num = 10, r = 0.15, deltaT = .000001, deltaX = .0025, nframes = 20, Q = q) 
  sampsNQ <- ifelse(sampsNQ < 0, 0, sampsNQ)
  sampsNQ <- ifelse(sampsNQ > max(bigD[,,j]), max(bigD[,,j]), sampsNQ)
  
  
  for(i in 1 : 10){
    step = begin + i
    compNQ[i,1] = i + begin                                      # the image i am attempting to predict
    compNQ[i,2] = ERROR(sampsNQ[,,i], bigD[,,(step)])            # Model error
    compNQ[i,3] = ERROR(bigD[,,begin], bigD[,,(step)])           # Persistance error
    compNQ[i,4] = compNQ[i,2] - compNQ[i,3]                      # Difference in errors
    compNQ[i,5] = (1 - (compNQ[i,2]/compNQ[i,3])) * 100          # Percent Change
  }
  write.csv(compNQ, file = paste("start_",begin,".csv", sep = ""), row.names = FALSE)
}




 # Create Q matrix and run the nowcast on begin timestep i
begin = 3
Q <- (bigD[,,begin] - bigD[,,(begin-1)])/100
sampsNQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                   timestep = begin, num = 10, r = 0.24, nframes = 20, Q = Q) 



for(i in 1 : 10){
  timestep = begin + i
  compNQ[i,1] = i + begin                                      # the image i am attempting to predict
  compNQ[i,2] = ERROR(sampsNQ[,,i], bigD[,,(timestep)])        # Model error
  compNQ[i,3] = ERROR(bigD[,,begin], bigD[,,(timestep)])       # Persistance error
  compNQ[i,4] = compNQ[i,2] - compNQ[i,3]                      # Difference in errors
  compNQ[i,5] = (1 - (compNQ[i,2]/compNQ[i,3])) * 100          # Percent Change
}




##########################################################################
### Normal Q data
##########################################################################
files = list.files('/Users/joshuanorth/Desktop/sfw/mses/', pattern = '*.csv',full.names = T)
Q = matrix(nrow = 60, ncol = 10)
for(i in 1: length(files)){
  temp = read.csv(files[i])
  temp2 = temp$Change
  Q[i,] = temp2
}

boxplot(Q, outline = F, xlab = "Timestep Projected", ylab = "Percent Increase in Accuracy",
        main = "Q")
abline(0,0)
abline(10,0)
abline(20,0)


files = list.files('/Users/joshuanorth/Desktop/mse_pos/', pattern = '*.csv',full.names = T)
Q = matrix(nrow = 60, ncol = 10)
for(i in 1: length(files)){
  temp = read.csv(files[i])
  temp2 = temp$Change
  Q[i,] = temp2
}

boxplot(Q, outline = F, xlab = "Timestep Projected", ylab = "Percent Increase in Accuracy",
        main = "Q")
abline(0,0)
abline(10,0)
abline(20,0)







files = list.files("/Users/joshuanorth/Desktop/MSE_q200_r01/", pattern = '*.csv',full.names = T)
Q200 = matrix(nrow = 60, ncol = 10)
for(i in 1: length(files)){
  temp = read.csv(files[i])
  temp2 = temp$Change
  Q200[i,] = temp2
}


boxplot(Q200, outline = F, xlab = "Timestep Projected", ylab = "Percent Increase in Accuracy",
        main = "Q69")
abline(0,0)
abline(10,0)
abline(20,0)



test = matrix(nrow = 56, ncol = 11)
test[1:56,1:10] = NoQ
NoQ = test

NoQ = NoQ[-c(14,15,16),]
Q = Q[-c(14,15,16),]
Qpos = Qpos[-c(14,15,16),]
Q200 = Q200[-c(14,15,16),]


png("Q_comparisons.png", height = 800, width = 1600)
boxplot(NoQ, at = c(1,6,11,16,21,26,31,36,41,46,51), outline = F, col = "purple", axes = F)
boxplot(Q, at =  c(2,7,12,17,22,27,32,37,42,47), outline = F, add = T, col = "green",
        xlab = "Timestep Projected", ylab = "Percent Increase in Accuracy",
        cex.lab = 1.5, cex.axis = 1.5)
boxplot(Qpos, at = c(3,8,13,18,23,28,33,38,43,48), outline = F, add = T, col = "gray", axes = F)
boxplot(Q200, at = c(4,9,14,19,24,29,34,39,44,49), outline = F, add = T, col = "red", axes = F)
abline(0,0)
legend("bottomright", title="Type of Q Used",
       c("None","Scaled by 1000","Scaled by 1000, Positive Only",
         "Scaled by 200"), fill=c("purple","green","gray","red"), horiz=F, cex=1.5)
dev.off()




