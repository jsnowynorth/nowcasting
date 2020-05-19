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

# get nowcasts ----------------------------------------------------------------

# timestep 10
begin = 10
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
est <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
               number = 10, timestep = begin, num = 10, 
               nframes = 20, Q = q, r = 0.15)
est10 <- ifelse(est < 0.0, 0, est)
est10 <- ifelse(est10 > 1.6, 1.6, est10)


# timestep 20
begin = 20
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
est <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
              number = 10, timestep = begin, num = 10, 
              nframes = 20, Q = q, r = 0.15)
est20 <- ifelse(est < 0.0, 0, est)
est20 <- ifelse(est20 > 1.6, 1.6, est20)


# timestep 30
begin = 30
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
est <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
               number = 10, timestep = begin, num = 10, 
               nframes = 20, Q = q, r = 0.15)
est30 <- ifelse(est < 0.0, 0, est)
est30 <- ifelse(est30 > 1.6, 1.6, est30)


# timestep 40
begin = 40
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
est <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
               number = 10, timestep = begin, num = 10, 
               nframes = 20, Q = q, r = 0.15)
est40 <- ifelse(est < 0.0, 0, est)
est40 <- ifelse(est40 > 1.6, 1.6, est40)



# timestep 50
begin = 50
# With Q
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
est <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
               number = 10, timestep = begin, num = 10, 
               nframes = 20, Q = q, r = 0.15)
est50 <- ifelse(est < 0.0, 0, est)
est50 <- ifelse(est50 > 1.6, 1.6, est50)




# plot nowcasts one -----------------------------------------------------------


mid_col = 'gray97'
png("/Users/joshuanorth/Desktop/f_vs_F.png", height = 2500, width = 1250)

layout.matrix <- matrix(c(1:10, 11, 11), nrow = 6, ncol = 2, byrow = T)

layout(mat = layout.matrix,
       heights = c(rep(1.5, 5), 0.5), # Heights of the two rows
       widths = rep(1.5, 2)) # Widths of the two columns

# layout.show(11)

par(mar = c(4,8,6,4),
    oma = c(4,1,1,1))


# first row
# 10
# par(mfrow = c(1,2))
image(est10[,,1] - bigD[,,10], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), main = "F(s,t) - f(s,t-1)", ylab = 't = 11',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,11] - bigD[,,10], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), main = "f(s,t) - f(s,t-1)",
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)

# par(mfrow = c(1,2))
# image.plot(est10[,,1] - bigD[,,10], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))
# image.plot(bigD[,,10] - bigD[,,9], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))


# second row
# 20
# par(mfrow = c(1,2))
image(est20[,,1] - bigD[,,20], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 21',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,21] - bigD[,,20], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)

# par(mfrow = c(1,2))
# image.plot(est20[,,1] - bigD[,,20], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))
# image.plot(bigD[,,20] - bigD[,,19], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))


# third row
# 30
# par(mfrow = c(1,2))
image(est30[,,1] - bigD[,,30], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 31',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,31] - bigD[,,30], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)

# par(mfrow = c(1,2))
# image.plot(est30[,,1] - bigD[,,30], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))
# image.plot(bigD[,,30] - bigD[,,29], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))


# fourth row
# 40
# par(mfrow = c(1,2))
image(est40[,,1] - bigD[,,40], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 41',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,41] - bigD[,,40], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)

# par(mfrow = c(1,2))
# image.plot(est40[,,1] - bigD[,,40], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))
# image.plot(bigD[,,40] - bigD[,,39], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))


# fifth row
# 50
# par(mfrow = c(1,2))
image(est50[,,1] - bigD[,,50], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 51',
           axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,51] - bigD[,,50], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
           axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)

# par(mfrow = c(1,2))
# image.plot(est50[,,1] - bigD[,,50], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))
# image.plot(bigD[,,50] - bigD[,,49], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue'))

# sixth row
image(x = seq(-1.6, 1.6, length.out = 256), z = t(t(seq(-1.6, 1.6, length.out = 256))), col = two.colors(start = 'blue', mid = mid_col), xlab='', axes = F)
axis(1, at = seq(-1.6, 1.6, 0.4), cex.axis = 4, outer = T, col = NA)

dev.off()


# plot nowcasts two -----------------------------------------------------------


mid_col = 'gray97'
png("/Users/joshuanorth/Desktop/f_vs_F_two.png", height = 2500, width = 1250)

layout.matrix <- matrix(c(1:10, 11, 11), nrow = 6, ncol = 2, byrow = T)

layout(mat = layout.matrix,
       heights = c(rep(1.5, 5), 0.5), # Heights of the two rows
       widths = rep(1.5, 2)) # Widths of the two columns

# layout.show(11)

par(mar = c(4,8,6,4),
    oma = c(4,1,1,1))


# first row
# 10
# par(mfrow = c(1,2))
image(est10[,,1] - bigD[,,10], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), main = "F(s,t+1) - f(s,t)", ylab = 't = 10',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,10] - bigD[,,9], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), main = "f(s,t) - f(s,t-1)",
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# second row
# 20
# par(mfrow = c(1,2))
image(est20[,,1] - bigD[,,20], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 20',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,20] - bigD[,,19], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# third row
# 30
# par(mfrow = c(1,2))
image(est30[,,1] - bigD[,,30], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 30',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,30] - bigD[,,29], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# fourth row
# 40
# par(mfrow = c(1,2))
image(est40[,,1] - bigD[,,40], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 40',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,40] - bigD[,,39], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# fifth row
# 50
# par(mfrow = c(1,2))
image(est50[,,1] - bigD[,,50], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 50',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,50] - bigD[,,49], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)

# sixth row
image(x = seq(-1.6, 1.6, length.out = 256), z = t(t(seq(-1.6, 1.6, length.out = 256))), col = two.colors(start = 'blue', mid = mid_col), xlab='', axes = F)
axis(1, at = seq(-1.6, 1.6, 0.4), cex.axis = 4, outer = T, col = NA)

dev.off()






















# plot nowcasts three -----------------------------------------------------------


mid_col = 'gray97'
png("/Users/joshuanorth/Desktop/f_vs_F_three.png", height = 2500, width = 1250)

layout.matrix <- matrix(c(1:10, 11, 11), nrow = 6, ncol = 2, byrow = T)

layout(mat = layout.matrix,
       heights = c(rep(1.5, 5), 0.5), # Heights of the two rows
       widths = rep(1.5, 2)) # Widths of the two columns

# layout.show(11)

par(mar = c(4,8,6,4),
    oma = c(4,1,1,1))


# first row
# 10
# par(mfrow = c(1,2))
image(bigD[,,11] - est10[,,1], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), main = "f(s,t) - F(s,t)", ylab = 't = 11',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,11] - bigD[,,10], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), main = "f(s,t) - f(s,t-1)",
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# second row
# 20
# par(mfrow = c(1,2))
image(bigD[,,21] - est20[,,1], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 21',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,21] - bigD[,,20], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# third row
# 30
# par(mfrow = c(1,2))
image(bigD[,,31] - est30[,,1], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 31',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,31] - bigD[,,30], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# fourth row
# 40
# par(mfrow = c(1,2))
image(bigD[,,41] - est40[,,1], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 41',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,41] - bigD[,,40], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)


# fifth row
# 50
# par(mfrow = c(1,2))
image(bigD[,,51] - est50[,,1], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col), ylab = 't = 51',
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)
image(bigD[,,51] - bigD[,,50], zlim = c(-1.6, 1.6), col = two.colors(start = 'blue', mid = mid_col),
      axes = F, cex.lab = 5, cex.main = 5, font.lab = 2)

# sixth row
image(x = seq(-1.6, 1.6, length.out = 256), z = t(t(seq(-1.6, 1.6, length.out = 256))), col = two.colors(start = 'blue', mid = mid_col), xlab='', axes = F)
axis(1, at = seq(-1.6, 1.6, 0.4), cex.axis = 4, outer = T, col = NA)

dev.off()




















