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

bigD <- load.netCDF.data(filepath = '~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/lightning_data')

run <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/Movement_Data.csv') # mac
run <- run[,2:5] # drop first column

move.perc <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/predict_data.csv', header = TRUE)
col.prec <- read.csv('~/Dropbox/Nowcasting/FDAWLM_USE_THIS_ONE/Final/data/column_predict.csv', header = TRUE)


head(run)


# vertical movement -------------------------------------------------------

plot(Horizontal ~ order, data = run, type = "b", xlab = "Time", ylab = "Number of Pizels Moved",
     main = "Horizontal Movement")


size.axis = 16
size.title = 18


ggplot(data = run, aes(y = Horizontal, x = order)) + 
  geom_point() +
  geom_line(linetype = "dotted") +
  geom_smooth(method='lm', colour = "black", se = FALSE) +
  xlab("Time") + 
  ylab("Number of Pixel's Moved") +
  ggtitle("Horizontal Movement") +
  theme(plot.title = element_text(size = size.title)) +
  theme(axis.title.x = element_text(size = size.axis)) + 
  theme(axis.title.y = element_text(size = size.axis))

ggplot(data = run, aes(y = Vertical, x = order)) + 
  geom_point() +
  geom_line(linetype = "dotted") +
  geom_smooth(method='lm', colour = "black", se = FALSE) +
  xlab("Time") + 
  ylab("Number of Pixel's Moved") +
  ggtitle("Vertical Movement") +
  theme(plot.title = element_text(size = size.title)) +
  theme(axis.title.x = element_text(size = size.axis)) + 
  theme(axis.title.y = element_text(size = size.axis))











