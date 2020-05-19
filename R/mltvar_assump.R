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
library(MVLM)

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


predict.regression <- function(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                               timestep = timestep, prev = prev, step = 10, run = run){
  
  # run <- points.move.partial(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = timestep)
  run <- run[(1:timestep),]
  a <- run[,1] # horizontal movement
  b <- run[,2] # vertical movement
  c <- run[,3] # error
  d <- run[,4] # order

  hor <- lm(a~d)
  ver <- lm(b~d)
  
  rho <- cor(run$Horizontal, run$Vertical)
  mod <- lm(cbind(Horizontal, Vertical) ~ order, data = run)
  mu1 <- t(matrix(c(1, 51))) %*% mod$coefficients[,1]
  mu2 <- t(matrix(c(1, 51))) %*% mod$coefficients[,2]
  
  
  
  hres <- mod$residuals[,1]
  vres <- mod$residuals[,2]
  plot((hres/sd(hres)), (vres/sd(vres)))
  
  
  if(is.matrix(prev) == TRUE){
    curr <- prev
  }else if(is.matrix(prev) == FALSE){
    curr <- bigD[,,timestep]
  }else{
    warning('Must input a valid timestep and or previous matrix')
  }
  
  tmp <- dim(curr)
  nx <- tmp[1]
  ny <- tmp[2]
  
  store <- array(NA, dim=c(400,400,step))
  
  for (k in 1:step){
    
    preds <- predict(mod, newdata = data.frame(order = (timestep+k)))
    horizontal = preds[,1]
    vertical = preds[,2]
    
    horizontal = predict(hor, newdata = data.frame(d = (timestep+k)))
    vertical = predict(ver, newdata = data.frame(d = (timestep+k)))
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
      
      nex[(ny+r_ver):ny,(nx-r_hor):nx] <- curr[(ny+r_ver):ny,(nx-r_hor):nx] 			# interpolate
      tmp <- m.right(obs1 = nex, obs2 = curr,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
      tmp_2 <- m.minus(obs1 = tmp$mright, obs2 = curr,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
      store[,,k] <- tmp_2$mminus 		# put into store matrix
      
    }else if(horizontal < 0 & vertical > 0){
      
      for(i in (1 - r_hor) : (nx)){
        for(j in (1 + r_ver) : ny){
          
          nex[j,i] <- curr[(j-r_ver),(i+r_hor)]
          
        } # end j for
      } # end i for
      
      nex[1 : r_ver, 1:-r_hor] <- curr[1 : r_ver, 1:-r_hor] 			# interpolate
      tmp <- m.left(obs1 = nex, obs2 = curr ,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
      tmp_2 <- m.plus(obs1 = tmp$mleft, obs2 = curr ,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
      store[,,k] <- tmp_2$mplus 		# put into store matrix
      
    }else{
      
      for(i in (1 - r_hor) : (nx)){
        for(j in 1 : (ny + r_hor) ){
          
          nex[j,i] <- curr[(j-r_ver),(i+r_hor)]
          
        } # end j for
      } # end i for
      
      nex[((ny + r_ver) : nx), 1 : -r_hor] <- curr[((ny + r_ver) : nx), 1 : -r_hor]  			# interpolate
      tmp <- m.left(obs1 = nex, obs2 = curr ,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
      tmp_2 <- m.minus(obs1 = tmp$mleft, obs2 = curr ,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
      store[,,k] <- tmp_2$mminus 		# put into store matrix
      
    }
  }
  
  return(store)
  
} # end function




# plot --------------------------------------------------------------------

# multivariate
begin = 50
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
mlt <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
              number = 10, timestep = begin, num = 10, 
              nframes = 20, Q = q, r = 0.15)


# univariate
univ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
               number = 10, timestep = begin, num = 10, 
               nframes = 20, Q = q, r = 0.15)



#Five time step staggered plot

mid_col = 'gray97'
# png("Q_Comparison_50", height = 2500, width = 1500)
png("/Users/joshuanorth/Desktop/mltv_Comparison_50.png", height = 2500, width = 2000)

layout.matrix <- matrix(c(1:20, 21, 21, 21, 22), nrow = 6, ncol = 4, byrow = T)

layout(mat = layout.matrix,
       heights = c(rep(1.5, 5), 0.5), # Heights of the two rows
       widths = rep(1.5, 4)) # Widths of the two columns

# layout.show(22)

par(mar = c(6,8,6,6),
    oma = c(4,1,1,1))
# first row
image(bigD[,,51], ylab = "51", main = "Actual", axes = F, col = tim.colors(256), cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(mlt[,,1], main = "Multivariate", axes = F, col = tim.colors(256), cex.main = 5, zlim = c(0, 1.6))
image(univ[,,1], main = "Univariate", axes = F, col = tim.colors(256), cex.main = 5, zlim = c(0, 1.6))
image(mlt[,,1] - univ[,,1], main = "Difference", axes = F, col = two.colors(start = 'blue', middle = mid_col), cex.main = 5, zlim = c(-1,1))

# second row
image(bigD[,,52], ylab = "52", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(mlt[,,2], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(univ[,,2], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(mlt[,,2] - univ[,,2], axes = F, col = two.colors(start = 'blue', middle = mid_col), zlim = c(-1,1))


# third row
image(bigD[,,53], ylab = "53", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(mlt[,,3], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(univ[,,3], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(mlt[,,3] - univ[,,3], axes = F, col = two.colors(start = 'blue', middle = mid_col), zlim = c(-1,1))

# fourth row
image(bigD[,,54], ylab = "54", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(mlt[,,4], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(univ[,,4], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(mlt[,,4] - univ[,,4], axes = F, col = two.colors(start = 'blue', middle = mid_col), zlim = c(-1,1))

# fifth row
image(bigD[,,55], ylab = "55", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(mlt[,,5], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(univ[,,5], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(mlt[,,5] - univ[,,5], axes = F, col = two.colors(start = 'blue', middle =mid_col), zlim = c(-1,1))

# sixth row
image(x = seq(0, 1.6, length.out = 256), z = t(t(seq(0, 1.6, length.out = 256))), col = tim.colors(256), xlab='', axes = F)
axis(1, at = seq(0, 1.6, 0.2), cex.axis = 4, outer = T, col = NA)

image(x = seq(-1, 1, length.out = 256), z = t(t(seq(-1, 1, length.out = 256))), col = two.colors(start = 'blue', middle = mid_col), xlab='', axes = F)
axis(1, at = seq(-1, 1, 0.5), cex.axis = 4, outer = T, col = NA)


dev.off()

