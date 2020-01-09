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
source('R/Load_netCDF_Data.R', chdir = TRUE)
source('R/Nowcasting_Function.R', chdir = TRUE)
source('R/Finite_Difference_Function.R', chdir = TRUE)
source('R/Move_Functions.R', chdir = TRUE)
source('R/Error_Function.R', chdir = TRUE)
source('R/Prediction_Functions.R', chdir = TRUE)


####################################################################
### Compile Data
####################################################################

bigD <- load.netCDF.data(filepath = 'data')

run <- read.csv('movement_data/Movement_Data.csv')
run <- run[,2:5] # drop first column

move.perc <- read.csv('movement_data/predict_data.csv', header = TRUE)
col.prec <- read.csv('movement_data/column_predict.csv', header = TRUE)



# How long move.all takes -------------------------------------------------

system.time(test <- move.all(bigD[,,50], bigD[,,51], tolerance = 0, iteration = 10, number = 10))


img <-  array(NA, dim = c(400,400,71))
x <- rep(NA, 71)
y <- rep(NA, 71)
er <- rep(NA,71)

system.time(
for(i in 1:71){
  tmp <- move.all(bigD[,,i], bigD[,,(i+1)], tolerance = 0, iteration = 10, number = 10)
  img[,,i] <- tmp$obs1
  x <- tmp$move[1]
  y <- tmp$move[2]
  er <- tmp$error
}
)



# Linear Regression -------------------------------------------------------

# horizontal
plot(run$Horizontal ~ run$order)

horModel = lm(Horizontal ~ order, data = run)
horBound = confint(horModel, level = 0.95)

x = seq(0,80, by = 1)
lowerHor = horBound[1] + horBound[2]*x
upperHor = horBound[3] + horBound[4]*x

abline(horModel)
points(lowerHor ~ x, type = 'l')
points(upperHor ~ x, type = 'l')

hor$coefficients <- c(horBound[1], horBound[2])

# vertical
plot(run$Vertical ~ run$order)

verModel = lm(Vertical ~ order, data = run)
verBound = confint(verModel, level = 0.95)

x = seq(0,80, by = 1)
lowerVer = verBound[1] + verBound[2]*x
upperVer = verBound[3] + verBound[4]*x

abline(verModel)
points(lowerVer ~ x, type = 'l')
points(upperVer ~ x, type = 'l')



# lower estimate
horBound <- c(horBound[1], horBound[2])
verBound <- c(verBound[1], verBound[2])

for(i in 35:68){
  
  # Initialize the begin timestep
  begin = i
  
  # Create an empty matrix to fill in with calculated values
  comp <- matrix(c(0), nc = 5, nr = 10)
  colnames(comp) <- c("Timestep", "FDAWLM", "Persisted", "Difference", "Change")
  
  # Create Q matrix and run the nowcast on begin timestep i
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  Q <- Q/200
  #Q = matrix(nrow = 400, ncol = 400, 0)
  samps <- nowcast.error(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = 50, num = 10, deltaT = .000001, deltaX = .0025, nframes = 25, hor = horBound, ver = verBound)
  
  for(i in 1 : 10){
    timestep = begin + i
    comp[i,1] = i + begin                                  # the image i am attempting to predict
    comp[i,2] = ERROR(samps[,,i], bigD[,,(timestep)])      # Model error
    comp[i,3] = ERROR(bigD[,,begin], bigD[,,(timestep)])   # Persistance error
    comp[i,4] = comp[i,2] - comp[i,3]                      # Difference in errors
    comp[i,5] = (1 - (comp[i,2]/comp[i,3])) * 100          # Percent Change
  }
  write.csv(comp, file = paste("start_",begin,".csv", sep = ""), row.names = FALSE)
}

# upper estimate
horBound <- c(horBound[3], horBound[4])
verBound <- c(verBound[3], verBound[4])


for(i in 21:22){
  
  # Initialize the begin timestep
  begin = i
  
  # Create an empty matrix to fill in with calculated values
  comp <- matrix(c(0), nc = 5, nr = 10)
  colnames(comp) <- c("Timestep", "FDAWLM", "Persisted", "Difference", "Change")
  
  # Create Q matrix and run the nowcast on begin timestep i
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  Q <- Q/200
  #Q = matrix(nrow = 400, ncol = 400, 0)
  samps <- nowcast.error(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = 50, num = 10, deltaT = .000001, deltaX = .0025, nframes = 25, hor = horBound, ver = verBound)
  
  for(i in 1 : 10){
    timestep = begin + i
    comp[i,1] = i + begin                                  # the image i am attempting to predict
    comp[i,2] = ERROR(samps[,,i], bigD[,,(timestep)])      # Model error
    comp[i,3] = ERROR(bigD[,,begin], bigD[,,(timestep)])   # Persistance error
    comp[i,4] = comp[i,2] - comp[i,3]                      # Difference in errors
    comp[i,5] = (1 - (comp[i,2]/comp[i,3])) * 100          # Percent Change
  }
  write.csv(comp, file = paste("start_",begin,".csv", sep = ""), row.names = FALSE)
}



# Improvement against persistance -----------------------------------------



system.time(trial <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = 10, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20))

comp <- matrix(c(0), nc = 4, nr = 10)
for(i in 1 : 10){
	comp[i,1] = ERROR(trial[,,i],bigD[,,(timestep+i)])
	comp[i,2] = ERROR(bigD[,,timestep],bigD[,,(timestep+i)])
	comp[i,3] = comp[i,1] - comp[i,2]
	comp[i,4] = 1 - (comp[i,2]/comp[i,1])
}
comp


# Rolling Window ----------------------------------------------------------


# timestep = n (start nowcast from)
# start = n (just a counter)
# write.csv(comp, file = "start_n", row.names = FALSE)

for(i in 41:58){
  
  # Initialize the begin timestep
  begin = i
  
  # Create an empty matrix to fill in with calculated values
  comp <- matrix(c(0), nc = 5, nr = 10)
  colnames(comp) <- c("Timestep", "FDAWLM", "Persisted", "Difference", "Change")
  
  # Create Q matrix and run the nowcast on begin timestep i
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  Q <- Q/200
  #Q = matrix(nrow = 400, ncol = 400, 0)
  samps <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = begin, num = 10, deltaT = .000001, deltaX = .0025, nframes = 25)
  
  for(i in 1 : 10){
    timestep = begin + i
    comp[i,1] = i + begin                                  # the image i am attempting to predict
    comp[i,2] = ERROR(samps[,,i], bigD[,,(timestep)])      # Model error
    comp[i,3] = ERROR(bigD[,,begin], bigD[,,(timestep)])   # Persistance error
    comp[i,4] = comp[i,2] - comp[i,3]                      # Difference in errors
    comp[i,5] = (1 - (comp[i,2]/comp[i,3])) * 100          # Percent Change
  }
  write.csv(comp, file = paste("start_",begin,".csv", sep = ""), row.names = FALSE)
}

par(mfrow = c(2,2))
image.plot(bigD[,,24])
image.plot(samps[,,1])
image.plot(bigD[,,33])
image.plot(samps[,,10])

# Run tests ---------------------------------------------------------------

smpl1 <- move.all(obs1 = bigD[,,30], obs2 = bigD[,,31], tolerance = 0, iteration = 10, number = 10)
smpl2 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = smpl1$obs1, nframes = 20)

first <- predict.regression(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = 30, prev = 29, step = 1, run = run)
second <- finiteD(deltaT = .000001, deltaX = .0025, uInit = first[,,1], nframes = 20)

ERROR(bigD[,,30],bigD[,,31])
ERROR(second, bigD[,,31])

smpl3 <- move.all(obs1 = bigD[,,44], obs2 = bigD[,,45], tolerance = 0, iteration = 10, number = 10)
smpl4 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = smpl3$obs1, nframes = 10)

system.time(smpl5 <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10))
system.time(smpl6 <- finiteD(deltaT = .000001, deltaX = .0025, uInit = smpl5$obs1, nframes = 20))

ERROR(bigD[,,60],bigD[,,61])
ERROR(smpl6, bigD[,,61])

image.plot(smpl6)

### Plot results and look at error

par(mfrow = c(1,2))
image.plot(abs(smpl2 - bigD[,,31]), zlim = c(0,2), main = "FDAWLM")
image.plot(abs(bigD[,,30] - bigD[,,31]), zlim = c(0,2), main = "Persistence")

ERROR(smpl2, bigD[,,31])
ERROR(bigD[,,30],bigD[,,31])


par(mfrow = c(1,2))
image.plot(abs(smpl4 - bigD[,,45]), zlim = c(0,2), main = "FDAWLM")
image.plot(abs(bigD[,,44] - bigD[,,45]), zlim = c(0,2), main = "Persistence")

ERROR(smpl4,bigD[,,45])
ERROR(bigD[,,44],bigD[,,45])


par(mfrow = c(1,2))
image.plot(abs(smpl6 - bigD[,,61]), zlim = c(0,2), main = "FDAWLM")
image.plot(abs(bigD[,,60] - bigD[,,61]), zlim = c(0,2), main = "Persistence")

ERROR(smpl6,bigD[,,61])
ERROR(bigD[,,60],bigD[,,61])



# Experiment with Q -------------------------------------------------------

Q <- matrix(c(0), nc = 400, nr = 400)
t_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
t <- finiteD(deltaT = .000001, deltaX = .0025, uInit = t_a$obs1, nframes = 20)

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


# Error Table, Q matrix comparison ----------------------------------------

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





 
 

# Some Plots --------------------------------------------------------------


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
 
 

#  Plots of the percent change data ---------------------------------------

files = list.files('/Users/joshuanorth/Dropbox/Nowcasting/Nowcasting/data/Qpositive', pattern = '*.csv',full.names = T)
Qpos = matrix(nrow = 56, ncol = 10)
for(i in 1: length(files)){
  temp = read.csv(files[i])
  temp2 = temp$Change
  Qpos[i,] = temp2
}

files = list.files('/Users/joshuanorth/Dropbox/Nowcasting/Nowcasting/data/Q', pattern = '*.csv',full.names = T)
Q = matrix(nrow = 56, ncol = 10)
for(i in 1: length(files)){
  temp = read.csv(files[i])
  temp2 = temp$Change
  Q[i,] = temp2
}

files = list.files('/Users/joshuanorth/Dropbox/Nowcasting/Nowcasting/data/NoQ', pattern = '*.csv',full.names = T)
NoQ = matrix(nrow = 56, ncol = 10)
for(i in 1: length(files)){
  temp = read.csv(files[i])
  temp2 = temp$Change
  NoQ[i,] = temp2
}

files = list.files('/Users/joshuanorth/Dropbox/Nowcasting/Nowcasting/data/Q200', pattern = '*.csv',full.names = T)
Q200 = matrix(nrow = 56, ncol = 10)
for(i in 1: length(files)){
  temp = read.csv(files[i])
  temp2 = temp$Change
  Q200[i,] = temp2
}




png("Movements", height = 800, width = 1600)

par(mfrow = c(1,2))
plot(run$Horizontal, type = "o", lty = 5, xlab = "Weather Timestep", 
     ylab = "Number of Pixel Moved", main = "Horizontal Movement", cex.lab = 1.5, cex.axis = 1.5,
     cex.main = 1.5)
abline(lm(run$Horizontal ~ run$order))

plot(run$Vertical, type = "o", lty = 5, xlab = "Weather Timestep",
     ylab = "Number of Pixel Moved", main = "Vertical Movement", cex.lab = 1.5, cex.axis = 1.5,
     cex.main = 1.5)
abline(lm(run$Vertical ~ run$order))

dev.off()



# Normal Q data -----------------------------------------------------------

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



# Staggered plot of nowcast to actual for timesteps 51-54 -----------------

  # Initialize the begin timestep
  begin = 50
  # Create Q matrix and run the nowcast on begin timestep i
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  Q <- Q/200
  samps1 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = begin, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20) 

  start_50 = samps1
  
  begin = 51
  # Create Q matrix and run the nowcast on begin timestep i
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  Q <- Q/200
  samps2 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = begin, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  start_51 = samps2
  
  begin = 52
  # Create Q matrix and run the nowcast on begin timestep i
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  Q <- Q/200
  samps3 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = begin, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  start_52 = samps3
  
  begin = 53
  # Create Q matrix and run the nowcast on begin timestep i
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  Q <- Q/200
  samps4 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = begin, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  start_53 = samps4
  
  
# Ten time step staggered plot

# png("staggered", height = 5000, width = 2500)
# 
# par(mfrow = c(10,5))
# par(mar = c(1,1,1,1))
# # first row
# image(bigD[,,51], ylab = "51", main = "Actual", axes = F, col = tim.colors(256),
#       cex.lab = 4, cex.main = 5, zlim = c(-1,2))
# image(start_50[,,1], main = "Project Starting at 51", axes = F, col = tim.colors(256),
#       cex.main = 5, zlim = c(-1,2))
# plot.new()
# plot.new()
# plot.new()
# # second row
# image(bigD[,,52], ylab = "52", axes = F, col = tim.colors(256),
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,1], main = "Project Starting at 52", axes = F, col = tim.colors(256),
#       cex.main = 5,zlim = c(-1,2))
# plot.new()
# plot.new()
# # third row
# image(bigD[,,53], ylab = "53", axes = F, col = tim.colors(256),
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,1], main = "Project Starting at 53", 
#       cex.main = 5, axes = F, col = tim.colors(256), zlim = c(-1,2))
# plot.new()
# # fourth row
# image(bigD[,,54], ylab = "54", axes = F, col = tim.colors(256), 
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_53[,,1], main = "Project Starting at 54", axes = F, col = tim.colors(256),
#       cex.main = 5, zlim = c(-1,2))
# # fifth row
# image(bigD[,,55], ylab = "55", axes = F, col = tim.colors(256),
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,5], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_53[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
# # sixth row
# image(bigD[,,56], ylab = "56", axes = F, col = tim.colors(256),
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,6], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,5], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_53[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
# # seventh row
# image(bigD[,,57], ylab = "57", axes = F, col = tim.colors(256),
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,7], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,6], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,5], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_53[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))
# # eighth row
# image(bigD[,,58], ylab = "58", axes = F, col = tim.colors(256), 
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,8], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,7], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,6], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_53[,,5], axes = F, col = tim.colors(256), zlim = c(-1,2))
# # ninth row
# image(bigD[,,59], ylab = "59", axes = F, col = tim.colors(256), 
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,9], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,8], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,7], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_53[,,6], axes = F, col = tim.colors(256), zlim = c(-1,2))
# # tenth row
# image(bigD[,,60], ylab = "60", axes = F, col = tim.colors(256), 
#       cex.lab = 4, zlim = c(-1,2))
# image(start_50[,,10], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_51[,,9], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_52[,,8], axes = F, col = tim.colors(256), zlim = c(-1,2))
# image(start_53[,,7], axes = F, col = tim.colors(256), zlim = c(-1,2))
# 
# dev.off()

# five by five staggered plot  

png("staggered", height = 2500, width = 2500)
  
par(mfrow = c(5,5))
par(mar = c(6,6,6,6))
# first row
image(bigD[,,51], ylab = "Timestep 51", main = "Observed Weather", axes = F, col = tim.colors(256),
        cex.lab = 4, cex.main = 5, zlim = c(-1,2))
image(start_50[,,1], main = "Project 51-55", axes = F, col = tim.colors(256),
        cex.main = 5, zlim = c(-1,2))
plot.new()
plot.new()
plot.new()
# second row
image(bigD[,,52], ylab = "Timestep 52", axes = F, col = tim.colors(256),
        cex.lab = 4, zlim = c(-1,2))
image(start_50[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_51[,,1], main = "Project 52-55", axes = F, col = tim.colors(256),
        cex.main = 5,zlim = c(-1,2))
plot.new()
plot.new()
# third row
image(bigD[,,53], ylab = "Timestep 53", axes = F, col = tim.colors(256),
        cex.lab = 4, zlim = c(-1,2))
image(start_50[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_52[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_52[,,1], main = "Project 53-54", 
        cex.main = 5, axes = F, col = tim.colors(256), zlim = c(-1,2))
plot.new()
# fourth row
image(bigD[,,54], ylab = "Timestep 54", axes = F, col = tim.colors(256), 
        cex.lab = 4, zlim = c(-1,2))
image(start_50[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_51[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_52[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_53[,,1], main = "Project 54-55", axes = F, col = tim.colors(256),
        cex.main = 5, zlim = c(-1,2))
# fifth row
image(bigD[,,55], ylab = "Timestep 55", axes = F, col = tim.colors(256),
        cex.lab = 4, zlim = c(-1,2))
image(start_50[,,5], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_51[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_52[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(start_53[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
dev.off()
  



# Error analysis ----------------------------------------------------------

# write smoothed images to file
for(i in 3:68){
  current = bigD[,,i]
  Q <- bigD[,,i] - bigD[,,(i-1)]
  Q <- Q/200
  for(j in 1:10){
    current = finiteD(deltaT = .000001, deltaX = .0025, uInit = current, nframes = 20)
    write.csv(current, file = paste("file_", i, "_", j, ".csv", sep = ""))
    print(c(i,j))
  }
}

# access smoothed images

files = list.files('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/smoothed_data', pattern = '*.csv',full.names = T)
store <- array(NA, dim = c(400,400,66,10))
for(i in 1:60){
  for(j in 1:10){
    temp = read.csv(files[((10*i)+j)])
    temp = as.matrix(temp[,-1])
    store[,,i,j] = temp
  }
}


# file_i_j is image i smoothed j times (i.e. timesteps)



for(j in 3:63){

begin = j
errors <- matrix(0, nrow = 50, ncol = 10)
moved_mats <- array(NA, dim = c(400,400,50))

for(i in 1:50){
  
  mu = c((horModel$coefficients[1] + horModel$coefficients[2]*begin),
         (verModel$coefficients[1] + verModel$coefficients[2]*begin))
  sigma = matrix(c(1.257^2,0,0,0.2765^2), ncol = 2, nrow = 2)
  
  tmp <- rep(0, 10)
  for(k in 1:10){
    curr <- store[,,begin,k]
    samps <- predict.regression.error(curr, tolerance = 0, iteration = 10, 
                                      number = 10, begin, step = 1, mu, sigma)
    tmp[k] <- ERROR(samps$mat[,,1], bigD[,,(begin+k)])
  }
  errors[i,] <- tmp
}

write.csv(errors, file = paste("start_",begin,".csv", sep = ""), row.names = FALSE)

}


# determine the pixel error in the movements

begin = 50
tm <- matrix(NA, ncol = 160000, nrow = 10)

for(k in 1:10){
  moved_mats <- array(NA, dim = c(400,400,50))
  
  for(i in 1:50){
      mu = c((horModel$coefficients[1] + horModel$coefficients[2]*(begin+k)),
           (verModel$coefficients[1] + verModel$coefficients[2]*(begin+k)))
      sigma = matrix(c(1.257^2,0,0,0.2765^2), ncol = 2, nrow = 2)
      curr <- store[,,begin,k]
      samps <- predict.regression.error(curr, tolerance = 0, iteration = 10, 
                                        number = 10, begin, step = 1, mu, sigma)
      moved_mats[,,i]<-samps$mat[,,1]
  }

  obs <- array(c(bigD[,,begin+k]), dim = c(400, 400, 50))
  
  moved_use <- ifelse(moved_mats<=0, NA, moved_mats)
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[k,] <- apply(obs_use >= moved_use, c(1,2), mean)
}

png("60", height = 750, width = 1000)

par(mfrow = c(2,5))
for(i in 1:10){
  vals <- tm[i,which(tm[i,]>0, arr.ind = T)]
  vals <- vals[which(vals<1, arr.ind = T)]
  hist(vals, main = begin+i)
}

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
  Q <- Q/amt[i]
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                  timestep = 50, num = 1, deltaT = .000001, deltaX = .0025, nframes = 20) 
  
  moved_use <- ifelse(pred[,,1]<=0, NA, pred[,,1])
  obs_use <- ifelse(obs<=0, NA, obs)
  
  tm[i,] <- apply(obs_use >= moved_use, c(1,2), mean)
  er[i,2] <- ERROR(pred[,,1], bigD[,,51])
}

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q

# Q = 101 minimizes error above, grid search around 100 by 10's
amt <- seq(1, 200, by = 10)
tm <- matrix(NA, ncol = 160000, nrow = 20)

er <- matrix(NA, ncol = 2, nrow = 20)
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

amt[which(er[,2] == min(er[,2]), arr.ind = T)] # 101 minimizes Q


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







# Error analysis with Q and movement --------------------------------------

# Randomly move the projection, and for each projection, scale Q by a randomly generated
# number with mean set to the optimal value for the image and with sd 10 (for now)

# mean value to scale Q by:
# 10: 35
# 20: 64
# 30: 52
# 40: 35
# 50: 99
# avg = 57


# PIT Histograms ----------------------------------------------------------


# determine mu and sigma
horModel = lm(Horizontal ~ order, data = run)
verModel = lm(Vertical ~ order, data = run)

begin = 3  # choose starting time point
tm <- matrix(NA, ncol = (160000*50), nrow = 56)  # initialize empty matrix to store PIT values in

for(k in 3:58){
  moved_mats <- array(NA, dim = c(400,400,50))
  
  for(i in 1:50){
    Q <- (bigD[,,k]-bigD[,,(k-1)]) * rnorm(1, mean = 1/57, sd = 10)
    mu = c((horModel$coefficients[1] + horModel$coefficients[2]*(begin+k)),
           (verModel$coefficients[1] + verModel$coefficients[2]*(begin+k)))
    sigma = matrix(c(1.257^2,0,0,0.2765^2), ncol = 2, nrow = 2)
    curr <- bigD[,,begin + (k-3)]
    samps <- predict.regression.error(obs = curr, tolerance = 0, iteration = 10, 
                                      number = 50, begin, step = 10, mu, sigma)
    moved_mats[,,i]<-samps$mat[,,10]
  }
  
  obs <- array(c(bigD[,,k]), dim = c(400, 400, 50))
  
  moved_mats <- ifelse(moved_mats<=0, 0, moved_mats)
  
  count <- 1
  for(m in 1:400){
    for(n in 1:400){
      tmp <- obs[m,n,] == moved_mats[m,n,]
      indT <- which(tmp == TRUE, arr.ind = T)
      indF <- which(tmp == FALSE, arr.ind = T)

      tmp[indT] <- runif(length(indT), min = 1/51, max = 1)  #sample(number, length(indT), replace = T)
      tmp[indF] <- mean(obs[m,n,indF] >= moved_mats[m,n,indF])
      tm[(k-2),c(seq(from = count, to = count+49, by = 1))] <- tmp
      
      count <- count + 50
    }
  }

}

newTm <- c(tm[])
hist(newTm, main = "PIT of 50 minute forecasts", xlab = "Probability Integral Transformation")
hist(newTm, main = "PIT of 50 minute forecasts", xlab = "Probability Integral Transformation", freq = F)







# Timing for est. the one-step-ahead shift forecast -----------------------
time <- rep(NA, 67)
for(i in 1:71){
  tmp <- system.time(test <- move.all(bigD[,,i], bigD[,,(i+1)], tolerance = 0, iteration = 10, number = 10))
  time[i] <- tmp[3]
}




system.time(test <- move.all(bigD[,,43], bigD[,,44], tolerance = 0, iteration = 10, number = 10))

start <- 1
steps <- 69
img <-  array(NA, dim = c(400,400,(steps+1)))
x <- rep(NA, (steps+1))
y <- rep(NA, (steps+1))
er <- rep(NA,(steps+1))
j <- 1

system.time(
  for(i in start:(start+steps)){
    tmp <- move.all(bigD[,,i], bigD[,,(i+1)], tolerance = 0, iteration = 10, number = 10)
    # img[,,j] <- tmp$obs1
    # x[j] <- tmp$move[1]
    # y[j] <- tmp$move[2]
    # er[j] <- tmp$error
    j <- j+1
  }
)


# time for 1 step is user(4.674) system(1.508) elapsed(6.188)
# time for 2 step is user(13.419) system(4.501) elapsed(17.931)
     


# PIT Histograms ----------------------------------------------------------

# determine mu and sigma
horModel = lm(Horizontal ~ order, data = run)
verModel = lm(Vertical ~ order, data = run)

begin = 3  # choose starting time point
tm <- matrix(NA, ncol = (160000*50), nrow = 56)  # initialize empty matrix to store PIT values in

for(k in 3:58){
  moved_mats <- array(NA, dim = c(400,400,50))
  
  for(i in 1:50){
    Q <- (bigD[,,k]-bigD[,,(k-1)]) * rnorm(1, mean = 1/57, sd = 10)
    mu = c((horModel$coefficients[1] + horModel$coefficients[2]*(begin+k)),
           (verModel$coefficients[1] + verModel$coefficients[2]*(begin+k)))
    sigma = matrix(c(1.257^2,0,0,0.2765^2), ncol = 2, nrow = 2)
    curr <- bigD[,,begin + (k-3)]
    samps <- predict.regression.error(obs = curr, tolerance = 0, iteration = 10, 
                                      number = 50, begin, step = 10, mu, sigma)
    moved_mats[,,i]<-samps$mat[,,10]
  }
  
  obs <- array(c(bigD[,,k]), dim = c(400, 400, 50))
  
  moved_mats <- ifelse(moved_mats<=0, 0, moved_mats)
  
  count <- 1
  for(m in 1:400){
    for(n in 1:400){
      tmp <- obs[m,n,] == moved_mats[m,n,]
      indT <- which(tmp == TRUE, arr.ind = T)
      indF <- which(tmp == FALSE, arr.ind = T)
      
      tmp[indT] <- runif(length(indT), min = 1/51, max = 1)  #sample(number, length(indT), replace = T)
      tmp[indF] <- mean(obs[m,n,indF] >= moved_mats[m,n,indF])
      tm[(k-2),c(seq(from = count, to = count+49, by = 1))] <- tmp
      
      count <- count + 50
    }
  }
  
}

newTm <- c(tm[])
hist(newTm, main = "PIT of 50 minute forecasts", xlab = "Probability Integral Transformation")
hist(newTm, main = "PIT of 50 minute forecasts", xlab = "Probability Integral Transformation", freq = F)




# Construct prediction fields ---------------------------------------------

store.images <- list()

for(i in 59:62){
  begin = i
  # With Q
  q <- bigD[,,begin] - bigD[,,(begin-1)]
  q <- q/100
  pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                  number = 10, timestep = begin, num = 10, 
                  nframes = 20, Q = q, r = 0.15)
  field <- ifelse(pred < 0.0, 0, pred)
  field <- ifelse(field > 1.6, 1.6, field)
  store.images[[i]] <- field
}

forecasts <- store.images

save(forecasts, file = "forecasts.RData")






