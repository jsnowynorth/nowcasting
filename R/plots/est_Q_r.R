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






# Estimate Q and r --------------------------------------------------------

Qscale <- seq(0, 500, by = 50)
delR <- seq(0, 0.25, by = 0.05)
err <- expand.grid(Qscale, delR)
err$RMSE1 <- 0
err$RMSE2 <- 0
err$RMSE3 <- 0
err$RMSE4 <- 0
err$RMSE5 <- 0
err$RMSE6 <- 0
err$RMSE7 <- 0
err$RMSE8 <- 0
err$RMSE9 <- 0
err$RMSE10 <- 0
colnames(err) <- c("Qscale", "r", "RMSE1", "RMSE2", "RMSE3", "RMSE4", "RMSE5", "RMSE6", "RMSE7", "RMSE8", "RMSE9", "RMSE10")



for(i in 59:66){
  
  for(j in 3:58){
    
    begin <- j
    if(err[i,1] == 0){
      Q <- (bigD[,,begin] - bigD[,,(begin-1)])
      Q <- ifelse(Q < 0.01, 0, Q)
    }else{
      Q <- (bigD[,,begin] - bigD[,,(begin-1)])/err[i,1]
      Q <- ifelse(Q < 0.01, 0, Q)
    }
    
    proj <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10,
                    timestep = begin, num = 10, r = err[i,2], nframes = 20, Q = Q)
    
    tmp.err1 <- ERROR(proj[,,1], bigD[,,(begin+1)])
    tmp.err2 <- ERROR(proj[,,2], bigD[,,(begin+2)])
    tmp.err3 <- ERROR(proj[,,3], bigD[,,(begin+3)])
    tmp.err4 <- ERROR(proj[,,4], bigD[,,(begin+4)])
    tmp.err5 <- ERROR(proj[,,5], bigD[,,(begin+5)])
    tmp.err6 <- ERROR(proj[,,6], bigD[,,(begin+6)])
    tmp.err7 <- ERROR(proj[,,7], bigD[,,(begin+7)])
    tmp.err8 <- ERROR(proj[,,8], bigD[,,(begin+8)])
    tmp.err9 <- ERROR(proj[,,9], bigD[,,(begin+9)])
    tmp.err10 <- ERROR(proj[,,10], bigD[,,(begin+10)])
    err[i,3] <- tmp.err1 + err[i,3]
    err[i,4] <- tmp.err2 + err[i,4]
    err[i,5] <- tmp.err3 + err[i,5]
    err[i,6] <- tmp.err4 + err[i,6]
    err[i,7] <- tmp.err5 + err[i,7]
    err[i,8] <- tmp.err6 + err[i,8]
    err[i,9] <- tmp.err7 + err[i,9]
    err[i,10] <- tmp.err8 + err[i,10]
    err[i,11] <- tmp.err9 + err[i,11]
    err[i,12] <- tmp.err10 + err[i,12]
    print(c(i, j))
  }
  err[i,3] <- err[i,3]/53
  err[i,4] <- err[i,4]/53
  err[i,5] <- err[i,5]/53
  err[i,6] <- err[i,6]/53
  err[i,7] <- err[i,7]/53
  err[i,8] <- err[i,8]/53
  err[i,9] <- err[i,9]/53
  err[i,10] <- err[i,10]/53
  err[i,11] <- err[i,11]/53
  err[i,12] <- err[i,12]/53
  print(err[i,])
}


write.csv(err, "/Users/joshuanorth/Desktop/est_q_r.csv")

optim <- read.csv("/Users/joshuanorth/Desktop/est_q_r.csv")

vals <- apply(optim[,4:13], 1, mean)

optim[which.min(vals),]

optim.dat <- data.frame(mse = vals, Q = err$Var1, r = err$Var2)


ggplot(data = optim.dat, aes(x = Q, y = r)) +
  geom_tile(aes(fill = mse)) +
  scale_fill_viridis(option = "plasma", limits = c(0.15,0.25)) +
  xlab("Q Values") +
  ylab("r Values")

ggplot(data = optim.dat, aes(x = Q, y = r)) +
  geom_point(aes(size = mse)) +
  scale_size_continuous(limits = c(0.17,0.25), name = "MSE") +
  xlab("Q Values") +
  ylab("r Values")
  




for(i in 1:dim(err)[1]){
  
  for(j in 3:62){
    
    begin <- j
    if(err[i,1] == 0){
      Q <- (bigD[,,begin] - bigD[,,(begin-1)])
      Q <- ifelse(Q < 0.01, 0, Q)
    }else{
      Q <- (bigD[,,begin] - bigD[,,(begin-1)])/err[i,1]
      Q <- ifelse(Q < 0.01, 0, Q)
    }
    
    proj <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 3,
                       timestep = begin, num = 3, r = err[i,2], nframes = 20, Q = Q)
    tmp.err1 <- ERROR(proj[,,1], bigD[,,(begin+1)])
    tmp.err2 <- ERROR(proj[,,2], bigD[,,(begin+2)])
    tmp.err3 <- ERROR(proj[,,3], bigD[,,(begin+3)])
    err[i,3] <- tmp.err1 + err[i,3]
    err[i,4] <- tmp.err2 + err[i,4]
    err[i,5] <- tmp.err3 + err[i,5]
    print(c(i, j))
  }
  err[i,3] <- err[i,3]/59
  err[i,4] <- err[i,4]/59
  err[i,5] <- err[i,5]/59
  print(err[i,])
}

write.csv(err, file = "/Users/joshuanorth/Desktop/err_1step/err_3.csv")














