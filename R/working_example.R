##----------------------------------------------------
## Name: Joshua North
##
## Date: 1/9/2019
##
## Project: Nowcasting
##
## Objective:
##
## Notes:
##
##----------------------------------------------------


# Load Libraries ----------------------------------------------------------
library(fields)
library(ncdf4)
library(RNetCDF)
library(xtable)
library(forecast)
library(plot3D)
library(ggplot2)
library(MASS)
library(mvtnorm)


# Souce Functions ---------------------------------------------------------
source('R/Load_netCDF_Data.R', chdir = TRUE)
source('R/Nowcasting_Function.R', chdir = TRUE)
source('R/Finite_Difference_Function.R', chdir = TRUE)
source('R/Move_Functions.R', chdir = TRUE)
source('R/Error_Function.R', chdir = TRUE)
source('R/Prediction_Functions.R', chdir = TRUE)

run <- read.csv('movement_data/Movement_Data.csv')
run <- run[,2:5]



# Load Data ---------------------------------------------------------------
bigD <- load.netCDF.data(filepath = 'data')

begin = 50
q <- bigD[,,begin] - bigD[,,(begin-1)]
q <- q/100
system.time(pred <- nowcast(bigD = bigD, tolerance = 0, iteration = 10, 
                number = 10, timestep = begin, num = 10, 
                nframes = 20, Q = q, r = 0.15))

image.plot(pred[,,2])










