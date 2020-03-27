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
library(cowplot)
library(ggpubr)
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



head(run)




# plots -------------------------------------------------------------------

# Initialize the begin timestep
begin = 50
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
wQ50 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
              number = 10, timestep = begin, num = 10, 
              nframes = 20, Q = q, r = 0.15)
wQ50 <- ifelse(wQ50 < 0.0, 0, wQ50)
wQ50 <- ifelse(wQ50 > 1.6, 1.6, wQ50)


# Initialize the begin timestep
begin = 51
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
wQ51 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                number = 10, timestep = begin, num = 10, 
                nframes = 20, Q = q, r = 0.15)
wQ51 <- ifelse(wQ51 < 0.0, 0, wQ51)
wQ51 <- ifelse(wQ51 > 1.6, 1.6, wQ51)

# Initialize the begin timestep
begin = 52
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
wQ52 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                number = 10, timestep = begin, num = 10, 
                nframes = 20, Q = q, r = 0.15)
wQ52 <- ifelse(wQ52 < 0.0, 0, wQ52)
wQ52 <- ifelse(wQ52 > 1.6, 1.6, wQ52)

# Initialize the begin timestep
begin = 53
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
wQ53 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                number = 10, timestep = begin, num = 10, 
                nframes = 20, Q = q, r = 0.15)
wQ53 <- ifelse(wQ53 < 0.0, 0, wQ53)
wQ53 <- ifelse(wQ53 > 1.6, 1.6, wQ53)

# Initialize the begin timestep
begin = 54
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
wQ54 <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                number = 10, timestep = begin, num = 10, 
                nframes = 20, Q = q, r = 0.15)
wQ54 <- ifelse(wQ50 < 0.0, 0, wQ54)
wQ54 <- ifelse(wQ50 > 1.6, 1.6, wQ54)



#Five time step staggered plot


layout.matrix <- matrix(c(1:31, 32,32,32,33), nrow = 7, ncol = 5, byrow = T)

# png("Run_50", height = 3000, width = 2500)
png("/Users/joshuanorth/Desktop/Run_50.png", height = 3000, width = 2500)

# par(mfrow = c(6,5))
layout(mat = layout.matrix,
       heights = c(rep(1.5, 6), 0.5), # Heights of the two rows
       widths = rep(1.5, 5)) # Widths of the two columns

# par(mar = c(6,8,6,6))
par(mar = c(6,8,6,6),
    oma = c(4,1,1,1))
# first row
image(bigD[,,51], ylab = "51", main = "Actual", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0,1.6), font.lab = 2)
image(wQ50[,,1], main = "Forecast Run at 50", axes = F, col = tim.colors(256),
      cex.main = 5, zlim = c(0,1.6))
plot.new()
plot.new()
plot.new()


# second row
image(bigD[,,52], ylab = "52", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0,1.6), font.lab = 2)
image(wQ50[,,2], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ51[,,1], main = "Forecast Run at 51", axes = F, col = tim.colors(256),
      cex.main = 5, zlim = c(0,1.6))
plot.new()
plot.new()


# third row
image(bigD[,,53], ylab = "53", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0,1.6), font.lab = 2)
image(wQ50[,,3], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ51[,,2], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ52[,,1], main = "Forecast Run at 52", axes = F, col = tim.colors(256),
      cex.main = 5, zlim = c(0,1.6))
plot.new()

# fourth row
image(bigD[,,54], ylab = "54", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0,1.6), font.lab = 2)
image(wQ50[,,4], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ51[,,3], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ52[,,2], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ53[,,1], main = "Forecast Run at 53", axes = F, col = tim.colors(256),
      cex.main = 5, zlim = c(0,1.6))

# fifth row
image(bigD[,,55], ylab = "55", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0,1.6), font.lab = 2)
image(wQ50[,,5], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ51[,,4], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ52[,,3], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ53[,,2], axes = F, col = tim.colors(256), zlim = c(0,1.6))

# sixth row
image(bigD[,,56], ylab = "56", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0,1.6), font.lab = 2)
image(wQ50[,,6], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ51[,,5], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ52[,,4], axes = F, col = tim.colors(256), zlim = c(0,1.6))
image(wQ53[,,3], axes = F, col = tim.colors(256), zlim = c(0,1.6))

# seventh row
plot.new()
image(x = seq(0, 1.6, length.out = 256), z = t(t(seq(0, 1.6, length.out = 256))), col = tim.colors(256), xlab='', axes = F)
axis(1, at = seq(0, 1.6, 0.2), cex.axis = 7, outer = T, col = NA)
plot.new()


dev.off()


















