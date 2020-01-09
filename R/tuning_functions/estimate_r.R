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

run <- read.csv('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Movement_Data.csv') 
run <- run[,2:5] # drop first column




# Estimate r --------------------------------------------------------------
  
Q <- matrix(c(0), nc = 400, nr = 400)
delR <- seq(0,1,by=.025)

delR.val <- NA
MSE <- NA

error.sum <- data.frame(MSE, delR.val)
ind <- 1

t_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
for(i in 1:length(delR)){
    t <- finiteD(r = delR[i], uInit = t_a$obs1, nframes = 20, Q)
    error.sum[ind,1] <- ERROR(t, bigD[,,61])
    error.sum[ind,2] <- delR[i]
    ind = ind+1
    print(ind)
}

error.sum[which.min(error.sum$MSE),]

# min mse is 0.1111874 at 0.25

# now a grid around 0.25

Q <- matrix(c(0), nc = 400, nr = 400)
delR <- seq(0.2,0.3,by=.001)

delR.val <- NA
MSE <- NA

error.sum <- data.frame(MSE, delR.val)
ind <- 1

t_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
for(i in 1:length(delR)){
  t <- finiteD(r = delR[i], uInit = t_a$obs1, nframes = 20, Q)
  error.sum[ind,1] <- ERROR(t, bigD[,,61])
  error.sum[ind,2] <- delR[i]
  ind = ind+1
  print(ind)
}

error.sum[which.min(error.sum$MSE),]

# min mse is 0.1111778 at 0.238

# now a grid around 0.238

Q <- matrix(c(0), nc = 400, nr = 400)
delR <- seq(0.21,0.24,by=.0005)

delR.val <- NA
MSE <- NA

error.sum <- data.frame(MSE, delR.val)
ind <- 1

t_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
for(i in 1:length(delR)){
  t <- finiteD(r = delR[i], uInit = t_a$obs1, nframes = 20, Q)
  error.sum[ind,1] <- ERROR(t, bigD[,,61])
  error.sum[ind,2] <- delR[i]
  ind = ind+1
  print(ind)
}

error.sum[which.min(error.sum$MSE),]

# min mse is 0.1111778 at 0.2375


# now a grid around 0.2375

q <- matrix(c(0), nc = 400, nr = 400)
delR <- seq(0.236,0.238,by=.00001)

delR.val <- NA
MSE <- NA

error.sum <- data.frame(MSE, delR.val)
ind <- 1

t_a <- move.all(obs1 = bigD[,,60], obs2 = bigD[,,61], tolerance = 0, iteration = 10, number = 10)
for(i in 1:length(delR)){
  t <- finiteD(r = delR[i], uInit = t_a$obs1, nframes = 20, Q = q)
  error.sum[ind,1] <- ERROR(t, bigD[,,61])
  error.sum[ind,2] <- delR[i]
  ind = ind+1
  print(ind)
}

error.sum[which.min(error.sum$MSE),]

# min mse is 0.1111778 at 0.2375



# testing it on multiple time steps

len <- 51
q <- matrix(c(0), nc = 400, nr = 400)
delR <- seq(0.001,0.25, length.out = len)
niter <- seq(5,60,by=5)
step.5 = step.10 = step.15 = step.20 = step.25 = step.30 = step.35 <- rep(NA, len)
step.40 = step.45 = step.50 = step.55 = step.60 <- rep(NA, len)
delR.val <- delR
error.sum <- data.frame(delR.val = delR, step.5, step.10, step.15, step.20, 
                        step.25, step.30, step.35, step.40, step.45, step.50, step.55, step.60)
tab.ind <- 2

for(l in niter){
  ind <- 1
  
  Q <- (bigD[,,l] - bigD[,,(l-1)])/69
  
  for(i in 1:length(delR)){
    
    pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 1, 
                    timestep = l, num = 1, r = delR[i], nframes = 20, Q = q)
    error.sum[ind,tab.ind] <- ERROR(pred[,,1], bigD[,,(l+1)])
    ind = ind+1
    print(c(ind,tab.ind))
  }
  
  tab.ind <- tab.ind + 1
}

error.sum2 <- apply(error.sum[,2:13], 1, mean)

plot(delR, error.sum2, type = "l",
     xlab = "Values of r", ylab = "MSE")

delR[which.min(error.sum2)]

# r = 0.03 is the best value












