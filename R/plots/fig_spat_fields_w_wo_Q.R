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

# plots -------------------------------------------------------------------

# Initialize the begin timestep
begin = 50
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
wQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                  number = 10, timestep = begin, num = 10, 
                  nframes = 20, Q = q, r = 0.15)
wQ <- ifelse(wQ < 0.0, 0, wQ)
wQ <- ifelse(wQ > 1.6, 1.6, wQ)


# Without Q
q <- matrix(rep(0, 160000), ncol = 400, nrow = 400)
nQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
              number = 10, timestep = begin, num = 10, 
              nframes = 20, Q = q, r = 0.15)
nQ <- ifelse(nQ < 0.01, 0, nQ)
nQ <- ifelse(nQ > 1.6, 1.6, nQ)



#Five time step staggered plot

mid_col = 'gray97'
# png("Q_Comparison_50", height = 2500, width = 1500)
png("/Users/joshuanorth/Desktop/Q_Comparison_50.png", height = 2500, width = 2000)

layout.matrix <- matrix(c(1:20, 21, 21, 21, 22), nrow = 6, ncol = 4, byrow = T)

layout(mat = layout.matrix,
       heights = c(rep(1.5, 5), 0.5), # Heights of the two rows
       widths = rep(1.5, 4)) # Widths of the two columns

# layout.show(22)

par(mar = c(6,8,6,6),
    oma = c(4,1,1,1))
# first row
image(bigD[,,51], ylab = "51", main = "Actual", axes = F, col = tim.colors(256), cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(wQ[,,1], main = "Projection with Q", axes = F, col = tim.colors(256), cex.main = 5, zlim = c(0, 1.6))
image(nQ[,,1], main = "Projection without Q", axes = F, col = tim.colors(256), cex.main = 5, zlim = c(0, 1.6))
image(wQ[,,1] - nQ[,,1], main = "Difference", axes = F, col = two.colors(start = 'blue', middle = mid_col), cex.main = 5, zlim = c(-1,1))

# second row
image(bigD[,,52], ylab = "52", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(wQ[,,2], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(nQ[,,2], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(wQ[,,2] - nQ[,,2], axes = F, col = two.colors(start = 'blue', middle = mid_col), zlim = c(-1,1))


# third row
image(bigD[,,53], ylab = "53", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(wQ[,,3], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(nQ[,,3], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(wQ[,,3] - nQ[,,3], axes = F, col = two.colors(start = 'blue', middle = mid_col), zlim = c(-1,1))

# fourth row
image(bigD[,,54], ylab = "54", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(wQ[,,4], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(nQ[,,4], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(wQ[,,4] - nQ[,,4], axes = F, col = two.colors(start = 'blue', middle = mid_col), zlim = c(-1,1))

# fifth row
image(bigD[,,55], ylab = "55", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(0, 1.6), font.lab = 2)
image(wQ[,,5], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(nQ[,,5], axes = F, col = tim.colors(256), zlim = c(0, 1.6))
image(wQ[,,5] - nQ[,,5], axes = F, col = two.colors(start = 'blue', middle =mid_col), zlim = c(-1,1))

# sixth row
image(x = seq(0, 1.6, length.out = 256), z = t(t(seq(0, 1.6, length.out = 256))), col = tim.colors(256), xlab='', axes = F)
axis(1, at = seq(0, 1.6, 0.2), cex.axis = 4, outer = T, col = NA)

image(x = seq(-1, 1, length.out = 256), z = t(t(seq(-1, 1, length.out = 256))), col = two.colors(start = 'blue', middle = mid_col), xlab='', axes = F)
axis(1, at = seq(-1, 1, 0.5), cex.axis = 4, outer = T, col = NA)


dev.off()



# plots -------------------------------------------------------------------

# Initialize the begin timestep
begin = 40
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
wQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
              number = 10, timestep = begin, num = 10, 
              nframes = 20, Q = q, r = 0.15)
wQ <- ifelse(wQ < 0.00, 0, wQ)
wQ <- ifelse(wQ > 1.6, 1.6, wQ)


# Without Q
q <- matrix(rep(0, 160000), ncol = 400, nrow = 400)
nQ <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
              number = 10, timestep = begin, num = 10, 
              nframes = 20, Q = q, r = 0.15)
nQ <- ifelse(nQ < 0.01, 0, nQ)
nQ <- ifelse(nQ > 1.6, 1.6, nQ)



#Five time step staggered plot

png("Q_Comparison_40", height = 2500, width = 1500)

par(mfrow = c(5,3))
par(mar = c(6,8,6,6))
# first row
image(bigD[,,41], ylab = "41", main = "Actual", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(-1,2), font.lab = 2)
image(wQ[,,1], main = "Projection with Q", axes = F, col = tim.colors(256),
      cex.main = 5, zlim = c(-1,2))
image(nQ[,,1], main = "Projection without Q", axes = F, col = tim.colors(256),
      cex.main = 5, zlim = c(-1,2))

# second row
image(bigD[,,42], ylab = "42", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(-1,2), font.lab = 2)
image(wQ[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(nQ[,,2], axes = F, col = tim.colors(256), zlim = c(-1,2))


# third row
image(bigD[,,43], ylab = "43", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(-1,2), font.lab = 2)
image(wQ[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(nQ[,,3], axes = F, col = tim.colors(256), zlim = c(-1,2))

# fourth row
image(bigD[,,44], ylab = "44", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(-1,2), font.lab = 2)
image(wQ[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(nQ[,,4], axes = F, col = tim.colors(256), zlim = c(-1,2))

# fifth row
image(bigD[,,45], ylab = "45", axes = F, col = tim.colors(256),
      cex.lab = 5, cex.main = 5, zlim = c(-1,2), font.lab = 2)
image(wQ[,,5], axes = F, col = tim.colors(256), zlim = c(-1,2))
image(nQ[,,5], axes = F, col = tim.colors(256), zlim = c(-1,2))


dev.off()
















