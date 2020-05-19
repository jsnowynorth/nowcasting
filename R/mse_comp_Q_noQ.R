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
# source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Load_netCDF_Data.R', chdir = TRUE)
# source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Nowcasting_Function.R', chdir = TRUE)
# source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Finite_Difference_Function.R', chdir = TRUE)
# source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Move_Functions.R', chdir = TRUE)
# source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Error_Function.R', chdir = TRUE)
# source('/Users/joshuanorth/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Prediction_Functions.R', chdir = TRUE)

source('R/Load_netCDF_Data.R', chdir = TRUE)
source('R/Nowcasting_Function.R', chdir = TRUE)
source('R/Finite_Difference_Function.R', chdir = TRUE)
source('R/Move_Functions.R', chdir = TRUE)
source('R/Error_Function.R', chdir = TRUE)
source('R/Prediction_Functions.R', chdir = TRUE)


####################################################################
### Compile Data
####################################################################

# bigD <- load.netCDF.data(filepath = '~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/lightning_data')
# 
# run <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Movement_Data.csv') # mac
# run <- run[,2:5] # drop first column
# 
# move.perc <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/predict_data.csv', header = TRUE)
# col.prec <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/column_predict.csv', header = TRUE)

bigD <- load.netCDF.data(filepath = 'data')

run <- read.csv('movement_data/Movement_Data.csv')
run <- run[,2:5] # drop first column

move.perc <- read.csv('movement_data/predict_data.csv', header = TRUE)
col.prec <- read.csv('movement_data/column_predict.csv', header = TRUE)

# compare Q values --------------------------------------------------------

compQ <- array(NA, dim = c(10, 4, 70))
colnames(compQ) <- c("FDAWLM", "Persisted", "Difference", "Change")


begin = 4
Q0 <- array(0, dim = dim(bigD[,,begin]))
proj <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10,
                timestep = (begin - 1), num = 1, r = 0, nframes = 5, Q = Q0)
q <- (bigD[,,begin] - proj[,,1])/50

for(j in 4:58){
  
  begin = j
  Q0 <- array(0, dim = dim(bigD[,,begin]))
  proj <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = 10,
                  timestep = (begin - 1), num = 1, r = 0.15, nframes = 5, Q = Q0)

  # q <- (bigD[,,begin] - bigD[,,(begin-1)])/150
  q <- (bigD[,,begin] - proj[,,1])/50
  sampQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                number = 10, timestep = begin, num = 10, 
                nframes = 20, Q = q, r = 0.15)
  
  
  for(i in 1 : 10){
    step = begin + i
    compQ[i,1,j] = ERROR(sampQ[,,i], bigD[,,(step)])^2              # Model error
    compQ[i,2,j] = ERROR(bigD[,,begin], bigD[,,(step)])^2           # Persistance error
    compQ[i,3,j] = compQ[i,1,j] - compQ[i,2,j]                    # Difference in errors
    compQ[i,4,j] = (1 - (compQ[i,1,j]/compQ[i,2,j])) * 100        # Percent Change
  }
  
  # q <- (bigD[,,(begin+1)] - sampQ[,,1])/50

  print(j)
  
}
saveRDS(compQ, file = "/Users/joshuanorth/Desktop/Qmse_reviewer.Rdata")




# No Q --------------------------------------------------------------------


compNQ <- array(NA, dim = c(10, 4, 70))
colnames(compNQ) <- c("FDAWLM", "Persisted", "Difference", "Change")

for(j in 4:58){
  
  begin = j
  Q0 <- array(0, dim = dim(bigD[,,begin]))
  sampQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                   number = 10, timestep = begin, num = 10, 
                   nframes = 20, Q = Q0, r = 0.2)
  
  for(i in 1 : 10){
    step = begin + i
    compNQ[i,1,j] = ERROR(sampQ[,,i], bigD[,,(step)])^2                # Model error
    compNQ[i,2,j] = ERROR(bigD[,,begin], bigD[,,(step)])^2             # Persistance error
    compNQ[i,3,j] = compNQ[i,1,j] - compNQ[i,2,j]                    # Difference in errors
    compNQ[i,4,j] = (1 - (compNQ[i,1,j]/compNQ[i,2,j])) * 100        # Percent Change
  }
  
}
saveRDS(compNQ, file = "/Users/joshuanorth/Desktop/NoQmse2.Rdata")


# compare results ---------------------------------------------------------


QOG <- readRDS("/Users/joshuanorth/Desktop/QmseOG.Rdata")
QOG_tab <- apply(QOG, c(1,2), mean, na.rm = T)

# QNew <- readRDS("/Users/joshuanorth/Desktop/Qmse_new.Rdata")
QNew <- readRDS("/Users/joshuanorth/Desktop/Qmse_reviewer.Rdata")
QNew_tab <- apply(QNew, c(1,2), mean, na.rm = T)

NoQ <- readRDS("/Users/joshuanorth/Desktop/NoQmse.Rdata")
NoQ_tab <- apply(NoQ, c(1,2), mean, na.rm = T)

Q1 <- QOG_tab[,4]
Q2 <- QNew_tab[,4]
Q3 <- NoQ_tab[,4]

xtable(rbind(Q1, Q2, Q3))


# cur <- apply(compQ, c(1,2), mean, na.rm = T) # original
# cur2 <- apply(compQ, c(1,2), mean, na.rm = T) # wrong new
# cur3 <- apply(compQ, c(1,2), mean, na.rm = T) # new
apply(compQ, c(1,2), mean, na.rm = T)
apply(compQ, c(1,2), quantile, probs = c(0.025, 0.975), na.rm = T)
apply(compNQ, c(1,2), mean, na.rm = T)
















