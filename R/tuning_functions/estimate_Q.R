####################################################################
### Joshua North
### 7/13/2017
### Finite Difference Approximation with Linear Movement (FDAWLM) script
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

#move.perc <- read.csv('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/predict_data.csv', header = TRUE)
#col.prec <- read.csv('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/column_predict.csv', header = TRUE)



####################################################################
### Titillate Q
####################################################################


Q <- matrix(c(0), nc = 400, nr = 400)
t_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t_a$obs1, nframes = 20)

t <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t_a$obs1, nframes = 20, Q)



Q <- (bigD[,,60])/100
t1_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t1 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t1_a$obs1, nframes = 20)

Q <- (bigD[,,60])/1000
t2_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t2 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t2_a $obs1, nframes = 20)

Q <- (sin((bigD[,,60])/1000)*cos((bigD[,,60])/1000))
t3_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t3 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t3_a $obs1, nframes = 20)

Q <- (bigD[,,60])/100000
t4_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t4 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t4_a $obs1, nframes = 20)

Q <- (bigD[,,60])/10
t5_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t5 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t5_a$obs1, nframes = 20)

Q <- bigD[,,60] - bigD[,,59]
Q <- (Q>0)/1000
t6_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t6 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t6_a$obs1, nframes = 20)

Q <- bigD[,,61] - bigD[,,60]
Q <- (Q>0)/1000
t7_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t7 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t7_a$obs1, nframes = 20)

### Q = (bigD[,,60])/100
par(mfrow = c(1,2))
image.plot(t, main = "No Q")
image.plot(t1, main = "Q = (bigD[,,60])/100")

par(mfrow = c(2,2))
image.plot(t, main = "No Q")
image.plot(t1, main = "Q = (bigD[,,60])/100")
image.plot(abs(t - bigD[,,61]), main = "No Q Error Forcast")
image.plot(abs(t1 - bigD[,,61]), main = "Q = (bigD[,,60])/100 Error Forcast")

### Q = (bigD[,,60])/1000
par(mfrow = c(1,2))
image.plot(t, main = "No Q")
image.plot(t2, main = "Q = (bigD[,,60])/1000")

par(mfrow = c(2,2))
image.plot(t, main = "No Q")
image.plot(t2, main = "Q = (bigD[,,60])/1000")
image.plot(abs(t - bigD[,,61]), main = "No Q Error Forcast")
image.plot(abs(t2 - bigD[,,61]), main = "Q = (bigD[,,60])/1000 Error Forcast")

### Q = sin((bigD[,,60])/1000)*cos((bigD[,,60])/1000)
par(mfrow = c(1,2))
image.plot(t, main = "No Q")
image.plot(t3, main = "Q = sin((bigD[,,60])/1000)*cos((bigD[,,60])/1000)")

par(mfrow = c(2,2))
image.plot(t, main = "No Q")
image.plot(t3, main = "Q = sin((bigD[,,60])/1000)*cos((bigD[,,60])/1000)")
image.plot(abs(t - bigD[,,61]), main = "No Q Error Forcast")
image.plot(abs(t3 - bigD[,,61]), main = "Q = sin((bigD[,,60])/1000)*cos((bigD[,,60])/1000)")

### Q = (bigD[,,60])/100000
par(mfrow = c(1,2))
image.plot(t, main = "No Q")
image.plot(t4, main = "Q = (bigD[,,60])/100000")

par(mfrow = c(2,2))
image.plot(t, main = "No Q")
image.plot(t4, main = "Q = (bigD[,,60])/100000")
image.plot(abs(t - bigD[,,61]), main = "No Q Error Forcast")
image.plot(abs(t4 - bigD[,,61]), main = "Q = (bigD[,,60])/100000 Error Forcast")

### Q = (bigD[,,60])/10
par(mfrow = c(1,2))
image.plot(t, main = "No Q")
image.plot(t5, main = "Q = (bigD[,,60])/10")

par(mfrow = c(2,2))
image.plot(t, main = "No Q")
image.plot(t5, main = "Q = (bigD[,,60])/10")
image.plot(abs(t - bigD[,,61]), main = "No Q Error Forcast")
image.plot(abs(t5 - bigD[,,61]), main = "Q = (bigD[,,60])/10 Error Forcast")

### Q = bigD[,,60] - bigD[,,59]
par(mfrow = c(1,2))
image.plot(t, main = "No Q")
image.plot(t6, main = "Q = (bigD[,,60] - bigD[,,59])/1000")

par(mfrow = c(2,2))
image.plot(t, main = "No Q")
image.plot(t6, main = "Q = bigD[,,60] - bigD[,,59]")
image.plot(abs(t - bigD[,,61]), main = "No Q Error Forcast")
image.plot(abs(t6 - bigD[,,61]), main = "(bigD[,,60] - bigD[,,59])/1000 Error Forcast")

### Q = bigD[,,61] - bigD[,,60]
par(mfrow = c(1,2))
image.plot(t, main = "No Q")
image.plot(t7, main = "Q = (bigD[,,61] - bigD[,,60])/1000")

par(mfrow = c(2,2))
image.plot(t, main = "No Q")
image.plot(t7, main = "Q = bigD[,,61] - bigD[,,60]")
image.plot(abs(t - bigD[,,61]), main = "No Q Error Forecast")
image.plot(abs(t7 - bigD[,,61]), main = "(bigD[,,61] - bigD[,,60])/1000 Error Forecast")

##########################################################################
### Error Table, Q matrix comparison
##########################################################################

###### Q non-zero

# starting timestep
begin = 51

# Create an empty matrix to fill in with calculated values
comp <- matrix(c(0), nc = 5, nr = 10)
colnames(comp) <- c("Timestep", "FDAWLM", "Persisted", "Difference", "Change")

# Create Q matrix and run the nowcast on begin timestep i
Q <- bigD[,,begin] - bigD[,,(begin-1)]
Q <- (Q>0)/1000
samps <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = begin, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20) 

for(i in 1 : 10){
  timestep = begin + i
  comp[i,1] = i + begin                                  # the image i am attempting to predict
  comp[i,2] = ERROR(samps[,,i], bigD[,,(timestep)])      # Model error
  comp[i,3] = ERROR(bigD[,,begin], bigD[,,(timestep)])   # Persistance error
  comp[i,4] = comp[i,2] - comp[i,3]                      # Difference in errors
  comp[i,5] = (1 - (comp[i,2]/comp[i,3])) * 100          # Percent Change
}

###### Q with zeros (no generation)
compNQ <- matrix(c(0), nc = 5, nr = 10)
colnames(compNQ) <- c("Timestep", "FDAWLM", "Persisted", "Difference", "Change")

# Create Q matrix and run the nowcast on begin timestep i
Q <- matrix(0, nrow = 400, ncol = 400)
sampsNQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = begin, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20) 

for(i in 1 : 10){
  timestep = begin + i
  compNQ[i,1] = i + begin                                  # the image i am attempting to predict
  compNQ[i,2] = ERROR(sampsNQ[,,i], bigD[,,(timestep)])      # Model error
  compNQ[i,3] = ERROR(bigD[,,begin], bigD[,,(timestep)])   # Persistance error
  compNQ[i,4] = compNQ[i,2] - compNQ[i,3]                      # Difference in errors
  compNQ[i,5] = (1 - (compNQ[i,2]/compNQ[i,3])) * 100          # Percent Change
}

xtable(comp, digits = c(0,0,3,3,3,6))
xtable(compNQ, digits = c(0,0,3,3,3,6))

comparison <- matrix(c(0), nc = 4, nr = 10)
colnames(comparison) <- c("Timestep", "With Q", "Without Q", "Difference")
comparison[,1] <- comp[,1]
comparison[,2] <- comp[,5]
comparison[,3] <- compNQ[,5]
comparison[,4] <- comp[,5] - compNQ[,5]
xtable(comparison, digits = c(0,0,6,6,6))

par(mfrow = c(2,2))
image.plot(bigD[,,50])
image.plot(bigD[,,51])
image.plot(Q)







##########################################################################
### Some Plots
##########################################################################


png("one.png", height = 400, width = 400) 
par(mfrow = c(1,2), mar = c(4,2,4,5))
image.plot(bigD[,,40], main = "Timestep 40", labels = F, tick = F)
image.plot(bigD[,,50], main = "Timestep 50", labels = F, tick = F)
dev.off()
png("two.png", height = 400, width = 400) 
par(mfrow = c(1,2), mar = c(4,2,4,5))
image.plot(bigD[,,60], main = "Timestep 60", labels = F, tick = F)
image.plot(bigD[,,70], main = "Timestep 70", labels = F, tick = F)
dev.off()


##########################################################################
### Normal Q data
##########################################################################
par(mfrow = c(1,1))

png("Q_1000.png", height = 500, width = 1000)
boxplot(Q, outline = F, xlab = "Timestep Projected", ylab = "Percent Increase in Accuracy",
        main = "Q scaled by 1000")
abline(0,0)
dev.off()

png("Q_pos.png", height = 500, width = 1000)
boxplot(Qpos, outline = F, xlab = "Timestep Projected", ylab = "Percent Increase in Accuracy",
        main = "Q scaled by 1000, only positive")
abline(0,0)
dev.off()


png("NoQ.png", height = 500, width = 1000)
boxplot(NoQ, outline = F, xlab = "Timestep Projected", ylab = "Percent Increase in Accuracy",
        main = "No Q used")
abline(0,0)
dev.off()

png("Q_200.png", height = 500, width = 1000)
boxplot(Q200, outline = F, xlab = "Timestep Projected", ylab = "Percent Reduction in MSE",
        main = "Q scaled by 200")
abline(0,0)
dev.off()


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


# Optimizing Q ------------------------------------------------------------

#### For timestep 50

# Grid search for best Q, by 100's
amt <- seq(1,2001, by = 100)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,50] - bigD[,,49]
  obs <- bigD[,,50]
  q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 50, num = 1, r = 0.03, nframes = 20, Q = q) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,51])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q

# Q = 101 minimizes error above, grid search around 100 by 10's
amt <- seq(90, 110, by = 1)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,50] - bigD[,,49]
  obs <- bigD[,,50]
  q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 50, num = 1, r = 0.2375, nframes = 20, Q = q) 
  
  # moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  # obs_use <- ifelse(obs<=0, NA, obs)
  
  moved_use <- ifelse(pred[,,1]=0, 0, pred[,,1])
  obs_use <- ifelse(obs=0, 0, obs)
  
  # moved_use <- pred[,,1]
  # obs_use <- obs
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,51])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q


pit <- rep(NA, 160000)
idx <- 1
for(i in 1:400){
  for(j in 1:400){
    if((obs_use[i,j] == 0) & (moved_use[i,j] == 0)){
      pit[idx] <- runif(1,0,1)
    }else{
      pit[idx] <- mean(obs_use[i,j] >= moved_use[i,j])
    }
    idx <- idx + 1
  }
}
hist(pit)



# Q = 101 minimizes error above, grid search around 100 by 1's
amt <- seq(90, 110, by = 1)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,50] - bigD[,,49]
  obs <- bigD[,,50]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 50, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,51])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)]

# this gives the optimal Q to be 99



#### For timestep 40

# Grid search for best Q, by 100's
amt <- seq(1,2001, by = 100)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,61] - bigD[,,60]
  obs <- bigD[,,61]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 60, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,62])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q


# Q = 101 minimizes error above, grid search around 100 by 10's
amt <- seq(1, 200, by = 10)
tm <- matrix(NA, ncol = 160000, nrow = 20)

er <- matrix(NA, ncol = 2, nrow = 20)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,40] - bigD[,,39]
  obs <- bigD[,,40]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 40, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,41])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 31 minimizes Q


# Q = 31 minimizes error above, grid search around 31 by 1's
amt <- seq(20, 40, by = 1)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,40] - bigD[,,39]
  obs <- bigD[,,40]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 40, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,41])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)]

# this gives the optimal Q to be 35



#### For timestep 30

# Grid search for best Q, by 100's
amt <- seq(1,2001, by = 100)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,30] - bigD[,,29]
  obs <- bigD[,,30]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 30, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,31])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q


# Q = 101 minimizes error above, grid search around 100 by 10's
amt <- seq(1, 200, by = 10)
tm <- matrix(NA, ncol = 160000, nrow = 20)

er <- matrix(NA, ncol = 2, nrow = 20)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,30] - bigD[,,29]
  obs <- bigD[,,30]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 30, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,31])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 51 minimizes Q


# Q = 51 minimizes error above, grid search around 51 by 1's
amt <- seq(40, 60, by = 1)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,30] - bigD[,,29]
  obs <- bigD[,,30]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 30, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,31])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)]

# this gives the optimal Q to be 52



#### For timestep 20

# Grid search for best Q, by 100's
amt <- seq(1,2001, by = 100)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,20] - bigD[,,19]
  obs <- bigD[,,20]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 20, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,21])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q


# Q = 101 minimizes error above, grid search around 100 by 10's
amt <- seq(1, 200, by = 10)
tm <- matrix(NA, ncol = 160000, nrow = 20)

er <- matrix(NA, ncol = 2, nrow = 20)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,20] - bigD[,,19]
  obs <- bigD[,,20]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 20, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,21])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 61 minimizes Q


# Q = 61 minimizes error above, grid search around 61 by 1's
amt <- seq(50, 70, by = 1)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,20] - bigD[,,19]
  obs <- bigD[,,20]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 20, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,21])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)]

# this gives the optimal Q to be 64



#### For timestep 10

# Grid search for best Q, by 100's
amt <- seq(1,2001, by = 100)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,10] - bigD[,,9]
  obs <- bigD[,,10]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 10, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,11])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q


# Q = 101 minimizes error above, grid search around 100 by 10's
amt <- seq(1, 200, by = 10)
tm <- matrix(NA, ncol = 160000, nrow = 20)

er <- matrix(NA, ncol = 2, nrow = 20)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,10] - bigD[,,9]
  obs <- bigD[,,10]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 10, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,11])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 31 minimizes Q


# Q = 31 minimizes error above, grid search around 31 by 1's
amt <- seq(20, 40, by = 1)
tm <- matrix(NA, ncol = 160000, nrow = 21)

er <- matrix(NA, ncol = 2, nrow = 21)
er[,1] <- amt

for(i in 1:length(amt)){
  Q <- bigD[,,10] - bigD[,,9]
  obs <- bigD[,,10]
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 10, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,11])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)]

# this gives the optimal Q to be 35




# Grid search Q -----------------------------------------------------------
len <- 51
Q.vals <- seq(1,201, length.out = len)
niter <- seq(5,60,by=5)
step.5 = step.10 = step.15 = step.20 = step.25 = step.30 = step.35 <- rep(NA, len)
step.40 = step.45 = step.50 = step.55 = step.60 <- rep(NA, len)
scaleQ <- delR
error.sum <- data.frame(delR.val = delR, step.5, step.10, step.15, step.20, 
                        step.25, step.30, step.35, step.40, step.45, step.50, step.55, step.60)

tab.ind <- 1

for(l in niter){
  Q <- bigD[,,l] - bigD[,,(l-1)]
  obs <- bigD[,,l]
  
  ind = 1
  tab.ind = tab.ind + 1
  
  for(i in 1:length(Q.vals)){
    q <- Q/Q.vals[i]
    pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                    timestep = l, num = 1, r = 0.03, nframes = 20, Q = q)
    
    error.sum[ind,tab.ind] <- ERROR(pred[,,1], bigD[,,(l+1)])
    
    print(c(ind,tab.ind))
    ind <- ind+1
    
  }
}


error.sum.q <- error.sum

error.sum2.q <- apply(error.sum.q[,2:13], 1, mean)

plot(Q.vals[5:51], error.sum2.q[5:51], type = "l", 
     xlab = "Scaling Values for Q", ylab = "MSE")

Q.vals[which.min(error.sum2.q)]

# Q scaled by 69 is best



