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
library(Matrix)

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

bigD <- load.netCDF.data(filepath = '~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/lightning_data')

run <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Movement_Data.csv') # mac
run <- run[,2:5] # drop first column

move.perc <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/predict_data.csv', header = TRUE)
col.prec <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/column_predict.csv', header = TRUE)



# edited move function ----------------------------------------------------



# movement function
predict.regression <- function(bigD = bigD, tolerance = 0, 
                               iteration = 10, number = 10, 
                               timestep = timestep, step = 10, 
                               hor.move, ver.move){

  
  curr <- bigD[,,timestep]
  curr <- bigD[,,50]
  
  tmp <- dim(curr)
  nx <- tmp[1]
  ny <- tmp[2]
  
  store <- array(NA, dim=c(400,400,step))
  
  for ( k in 1:step){
    
    horizontal = hor.move
    vertical = ver.move
    r_hor <- floor(horizontal)
    r_ver <- floor(vertical)
    hor_remainder <- horizontal - r_hor
    ver_remainder <- vertical - r_ver
    
    nex <- matrix(c(0), nr = nx, nc = ny)
    
    if(horizontal > 0 & vertical > 0){
      
      for(i in 1 : (nx - r_hor) ){			# choose 1 : nx-r_hor (whole number) of matrix
        for(j in (1 + r_ver) : ny ){		# choose 1+r_ver (whole number) of matrix
          
          nex[j,i] <- curr[(j-r_ver),(i+r_hor)]		# shift matrix down and right (moves image up and right)
          
        } # end for j
      } # end for i	
      
      nex[(1:r_hor),((ny - r_ver) : ny)] <- curr[(1:r_hor),((ny - r_ver) : ny)]			# interpolate
      tmp <- m.right(obs1 = nex, obs2 = curr,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
      tmp_2 <- m.plus(obs1 = tmp$mright, obs2 = curr,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
      store[,,k] <- tmp_2$mplus 		# put into store matrix
      
    }else if(horizontal > 0 & vertical < 0){
      
      for(i in 1 : (nx - r_hor) ){ 			# move integer value
        for(j in 1 : (ny + r_ver) ){
          
          nex[j,i] <- curr[(j-r_ver),(i+r_hor)]
          
        } # end for j
      } # end for i 
      
      nex[(1:r_hor),(1:-r_ver)] <- curr[(1:r_hor),(1:-r_ver)] 			# interpolate
      tmp <- m.right(obs1 = nex, obs2 = curr,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
      tmp_2 <- m.minus(obs1 = tmp$mright, obs2 = curr,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
      store[,,k] <- tmp_2$mminus 		# put into store matrix
      
    }else if(horizontal < 0 & vertical > 0){
      
      for(i in (1 - r_hor) : (nx + r_hor)){
        for(j in (1 + r_ver) : ny){
          
          nex[j,i] <- curr[(j-r_ver),(i-r_hor)]
          
        } # end j for
      } # end i for
      
      nex[((nx + r_hor) : nx),((ny - r_ver) : ny)] <- curr[((nx + r_hor) : nx),((ny - r_ver) : ny)]  			# interpolate
      tmp <- m.left(obs1 = nex, obs2 = curr ,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
      tmp_2 <- m.plus(obs1 = tmp$mleft, obs2 = curr ,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
      store[,,k] <- tmp_2$mplus 		# put into store matrix
      
    }else{
      
      for(i in (1 - r_hor) : (nx + r_hor)){
        for(j in (1 - r_ver) : (ny + r_ver) ){
          
          nex[j,i] <- curr[(j-r_ver),(i-r_hor)]
          
        } # end j for
      } # end i for
      
      nex[((nx + r_hor) : nx),((ny + r_ver) : ny)] <- curr[((nx + r_hor) : nx),((ny + r_ver) : ny)]  			# interpolate
      tmp <- m.left(obs1 = nex, obs2 = curr ,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
      tmp_2 <- m.minus(obs1 = tmp$mleft, obs2 = curr ,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
      store[,,k] <- tmp_2$mminus 		# put into store matrix
      
    }
  }
  
  return(store)
  
} # end function


nowcast <- function(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                    timestep = timestep, num = 10, r = 0.2375, deltaT = .000001, 
                    deltaX = .0025, nframes = 20, hor.move, ver.move, Q){
  
  store <- array(NA, dim=c(400,400,num))
  prev = 0
  
  for(i in 1 : num){
    
    smpl1 <- predict.regression(bigD = bigD, tolerance = 0, 
                                iteration = 10, number = 10, 
                                timestep = timestep, step = 10, 
                                hor.move, ver.move)
    temp <- smpl1[,,1]
    smpl2 <- finiteD(r, deltaT, deltaX, uInit = temp, nframes, Q)
    
    store[,,i] <- smpl2
    prev <- smpl2
    
  } # end for
  
  return(store = store)
  
} # end nowcast

frct1 <- matrix(0, ncol = 100, nrow = 160000)
frct2 <- matrix(0, ncol = 100, nrow = 160000)
frct3 <- matrix(0, ncol = 100, nrow = 160000)
frct4 <- matrix(0, ncol = 100, nrow = 160000)
frct5 <- matrix(0, ncol = 100, nrow = 160000)
begin = 50
ahd <- 5

for(i in 1:100){
  hor.mod <- lm(Horizontal ~ order, data = run)
  hor.step = predict(hor.mod, newdata = data.frame(order = (begin+1)))
  hm = rnorm(n = 1, mean = hor.step, sd = summary(hor.mod)$sigma)
  
  ver.mod <- lm(Vertical ~ order, data = run)
  ver.step = predict(ver.mod, newdata = data.frame(order = (begin+1)))
  vm = rnorm(n = 1, mean = ver.step, sd = summary(ver.mod)$sigma)
  
  Q <- bigD[,,begin] - bigD[,,(begin-1)]
  q <- Q/100
  trial = nowcast(bigD = bigD, tolerance = 0, iteration = 10, number = ahd, 
                  timestep = begin, num = ahd, r = 0.15, nframes = 20, hor.move = hm, ver.move = vm, Q = q)
  
  trial <- ifelse(trial<0.01, 0, trial)
  trial <- ifelse(trial>1.6, 1.6, trial)
  frct1[,i] = c(trial[,,1])
  frct2[,i] = c(trial[,,2])
  frct3[,i] = c(trial[,,3])
  frct4[,i] = c(trial[,,4])
  frct5[,i] = c(trial[,,5])
  print(i)
}


ary <- c(bigD[,,53])
rnk <- rep(NA, length(ary))
list = seq(1:100)

for(i in 1:length(ary)){
  tmp.num <- ary[i]
  tmp.len <- sum(frct3[i,] == tmp.num)
  if(tmp.num == 0){
    if(tmp.len == 0){
      rnk[i] <-1
    }else{
      rnk[i] <- sample(seq(1:(tmp.len+1)), size = 1)
    }
  }else if(tmp.num > 0){
    rnk[i] <- rank(c(ary[i], frct3[i,]), ties.method = "random")[1]
  }
  print(i)
}

hist((rnk/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F)




# PIT 51 ------------------------------------------------------------------

ary <- c(bigD[,,51])
rnk1 <- rep(NA, length(ary))
list = seq(1:100)

for(i in 1:length(ary)){
  tmp.num <- ary[i]
  tmp.len <- sum(frct1[i,] == tmp.num)
  if(tmp.num == 0){
    if(tmp.len == 0){
      rnk1[i] <-1
    }else{
      rnk1[i] <- sample(seq(1:(tmp.len+1)), size = 1)
    }
  }else if(tmp.num > 0){
    rnk1[i] <- rank(c(ary[i], frct1[i,]), ties.method = "random")[1]
  }
  print(i)
}

hist((rnk1/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F)



# PIT 52 ------------------------------------------------------------------

ary <- c(bigD[,,52])
rnk2 <- rep(NA, length(ary))
list = seq(1:100)

for(i in 1:length(ary)){
  tmp.num <- ary[i]
  tmp.len <- sum(frct2[i,] == tmp.num)
  if(tmp.num == 0){
    if(tmp.len == 0){
      rnk2[i] <-1
    }else{
      rnk2[i] <- sample(seq(1:(tmp.len+1)), size = 1)
    }
  }else if(tmp.num > 0){
    rnk2[i] <- rank(c(ary[i], frct2[i,]), ties.method = "random")[1]
  }
  print(i)
}

hist((rnk2/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F)




# PIT 53 ------------------------------------------------------------------

ary <- c(bigD[,,53])
rnk3 <- rep(NA, length(ary))
list = seq(1:100)

for(i in 1:length(ary)){
  tmp.num <- ary[i]
  tmp.len <- sum(frct3[i,] == tmp.num)
  if(tmp.num == 0){
    if(tmp.len == 0){
      rnk3[i] <-1
    }else{
      rnk3[i] <- sample(seq(1:(tmp.len+1)), size = 1)
    }
  }else if(tmp.num > 0){
    rnk3[i] <- rank(c(ary[i], frct3[i,]), ties.method = "random")[1]
  }
  print(i)
}

hist((rnk3/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F)




# PIT 54 ------------------------------------------------------------------

ary <- c(bigD[,,54])
rnk4 <- rep(NA, length(ary))
list = seq(1:100)

for(i in 1:length(ary)){
  tmp.num <- ary[i]
  tmp.len <- sum(frct4[i,] == tmp.num)
  if(tmp.num == 0){
    if(tmp.len == 0){
      rnk4[i] <-1
    }else{
      rnk4[i] <- sample(seq(1:(tmp.len+1)), size = 1)
    }
  }else if(tmp.num > 0){
    rnk4[i] <- rank(c(ary[i], frct4[i,]), ties.method = "random")[1]
  }
  print(i)
}

hist((rnk4/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F)




# PIT 55 ------------------------------------------------------------------

ary <- c(bigD[,,55])
rnk5 <- rep(NA, length(ary))
list = seq(1:100)

for(i in 1:length(ary)){
  tmp.num <- ary[i]
  tmp.len <- sum(frct5[i,] == tmp.num)
  if(tmp.num == 0){
    if(tmp.len == 0){
      rnk5[i] <-1
    }else{
      rnk5[i] <- sample(seq(1:(tmp.len+1)), size = 1)
    }
  }else if(tmp.num > 0){
    rnk5[i] <- rank(c(ary[i], frct5[i,]), ties.method = "random")[1]
  }
  print(i)
}

hist((rnk5/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F)


par(mfrow = c(1,5))
hist((rnk1/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F, ylim = c(0,5))
hist((rnk2/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F, ylim = c(0,5))
hist((rnk3/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F, ylim = c(0,5))
hist((rnk4/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F, ylim = c(0,5))
hist((rnk5/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F, ylim = c(0,5))


# Pit trials --------------------------------------------------------------


ranks = apply(obs, 1, rank, ties.method = "random")
hist(ranks)

ary <- c(bigD[,,50])
rnk <- rep(NA, length(ary))
list = c(1,2,3)
for(i in 1:length(ary)){
  if((ary[i] == 0 && obs[i,1:2] == 0)){
    rnk[i] <-  sample(c(1,2,3), size = 1)
  }else if(obs[i,1] > 0 && obs[i,2] == 0 && ary[i] == 0){
    rnk[i] <- sample(c(1,2), size = 1)
  }else if(obs[i,1:2] > 0 && ary[i] == 0){
    rnk[i] = 1
  }else{
    tmp.rnk <- rank(c(ary[i], obs[i,1:2]))
    rnk[i] <- tmp.rnk[1]
  }
  print(i)
}


ary <- c(bigD[,,50])
rnk <- NA
list = c(1,2,3)

ind <- which(which(obs[,1:2] == 0) %in% which(ary == 0))
length(ind)
rnk <- sample(list, replace = T, size = length(ind))

ind.2 <- which(which(obs[,1] != 0) %in% which(obs[,2] == 0) %in% which(ary == 0))
length(ind.2)
rnk <- c(rnk, sample(c(1,2), replace = T, size = length(ind.2)))

ind.3 <- which(which(obs[,2] != 0) %in% which(obs[,1] == 0) %in% which(ary == 0))
length(ind.3)
rnk <- c(rnk, sample(c(1,2), replace = T, size = length(ind.3)))

ind.4 <- which(which(obs[,1:2] != 0) %in% which(ary == 0))
length(ind.4)
rnk <- c(rnk, sample(c(1,2), replace = T, size = length(ind.4)))

ind.5 <- which(ary != 0)
length(ind.5)
rnk <- c(rnk, rep(1, length(ind.5)))
length(rnk)









ary <- c(bigD[,,50])
rnk <- NA
list = seq(1:100)

ind <- which(which(obs[,] == 0) %in% which(ary == 0))
length(ind)
rnk <- sample(list, replace = T, size = length(ind))

ind.4 <- which(which(obs[,] != 0) %in% which(ary == 0))
length(ind.4)
rnk <- c(rnk, rep(1, length(ind.4)))

ind.5 <- which(ary != 0)
length(ind.5)
rnks <- apply(cbind(ary[ind.5], obs[ind.5,]), 1, rank, ties.method = "random")
rnk <- c(rnk, rnks[1,])
length(rnk)

hist(rnk/100, main = "Probability Integral Transform, Time = 50", xlab = "Density", ylab = "Frequency", freq = F)





# attempt 2 ---------------------------------------------------------------

ary <- c(bigD[,,51])
rnk <- rep(NA, length(ary))
list = seq(1:100)

for(i in 1:length(ary)){
  tmp.num <- ary[i]
  tmp.len <- sum(obs[i,] == tmp.num)
  #tmp.len <- length(which(obs[i,] == tmp.num))
  if(tmp.num == 0){
    if(tmp.len == 0){
      rnk[i] <-1
    }else{
      rnk[i] <- sample(seq(1:(tmp.len+1)), size = 1)
    }
  }else if(tmp.num > 0){
    rnk[i] <- rank(c(ary[i], obs[i,]), ties.method = "random")[1]
  }
  print(i)
}

hist((rnk/101)-(1/202), main = "", xlab = "Density", ylab = "Randomized Rank", freq = F)



