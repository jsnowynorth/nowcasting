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
library(mvnormtest)
library(energy)
library(Hmisc)

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


# vertical movement -------------------------------------------------------

plot(Horizontal ~ order, data = run, type = "b", xlab = "Time", ylab = "Number of Pizels Moved",
     main = "Horizontal Movement")


size.axis = 16
size.title = 18


p1 <- ggplot(data = run, aes(y = Horizontal, x = order)) + 
  geom_point() +
  geom_line(linetype = "dotted") +
  geom_smooth(method='lm', colour = "black", se = FALSE) +
  xlab("Time") + 
  ylab("Number of Pixel's Moved") +
  ggtitle("Easting") +
  theme(plot.title = element_text(size = size.title)) +
  theme(axis.title.x = element_text(size = size.axis)) + 
  theme(axis.title.y = element_text(size = size.axis))

p2 <- ggplot(data = run, aes(y = Vertical, x = order)) + 
  geom_point() +
  geom_line(linetype = "dotted") +
  geom_smooth(method='lm', colour = "black", se = FALSE) +
  xlab("Time") + 
  ylab("Number of Pixel's Moved") +
  ggtitle("Northing") +
  theme(plot.title = element_text(size = size.title)) +
  theme(axis.title.x = element_text(size = size.axis)) + 
  theme(axis.title.y = element_text(size = size.axis))


ggsave("/Users/joshuanorth/Desktop/movement.png", cowplot::plot_grid(p1, p2, nrow = 1), width = 16, height = 6)

mod <- lm(cbind(Horizontal, Vertical) ~ order, data = run[-c(60:67),])
predict(run[c(60:67),], mod)
preds <- predict(mod, run[c(60:67),])
preds - run[c(60:67),c(1,2)]

h_res <- mod$residuals[,1]
v_res <- mod$residuals[,2]

h_mod_sd <- (h_res - mean(h_res)/sd(h_res))
v_mod_sd <- (v_res - mean(v_res)/sd(v_res))

plot(h_res ~ v_res, asp=1)
plot(h_mod_sd ~ v_mod_sd)
plot(h_res ~ v_res, asp = 1, xlab = 'Vertical Residuals', ylab = 'Horzontal Residuals', main = 'Standardized Residuals')
plot(h_res ~ v_res, xlab = 'Vertical Residuals', ylab = 'Horzontal Residuals', main = 'Standardized Residuals',
     xlim = c(-3, 3), ylim = c(-3,3))
abline(v = -2)
abline(v = 2)
abline(h = -2)
abline(h = 2)



h_mod <- lm(Horizontal ~ order, data = run)
v_mod <- lm(Vertical ~ order, data = run)

hres <- h_mod$residuals
vres <- v_mod$residuals

hsd <- hres/sd(hres)
vsd <- vres/sd(vres)

cor.test(hsd, vsd, method = 'pearson')

inds <- seq(1,67)[hsd < abs(2) & vsd < abs(2)]
hinds <- seq(1,67)[hsd < abs(2)]
vinds <- seq(1,67)[vsd < abs(2)]

plot(hsd, vsd, xlim = c(-3, 3), ylim = c(-3, 3))
plot(hsd[inds], vsd[inds], xlim = c(-3, 3), ylim = c(-3, 3))
cor.test(hsd[inds], vsd[inds], method = 'pearson')
cor.test(hsd[inds], vsd[inds], method = 'kendall')
cor.test(hsd[inds], vsd[inds], method = 'spearman')


plot(h_mod$residuals[inds], v_mod$residuals[inds])

cor.test(run$Horizontal, run$Vertical)
cor.test(h_mod$residuals, v_mod$residuals, method = 'pearson')
cor.test(h_mod$residuals[inds], v_mod$residuals[inds], method = 'pearson')
cor.test(h_mod$residuals, v_mod$residuals, method = 'kendall')
cor.test(h_mod$residuals, v_mod$residuals, method = 'spearman')








h_res <- rstandard(h_mod)
v_res <- rstandard(v_mod)

h_res <- rstudent(h_mod)
v_res <- rstudent(v_mod)

h_mod_sd <- (h_mod$residuals - mean(h_mod$residuals)/sd(h_mod$residuals))
v_mod_sd <- (v_mod$residuals - mean(v_mod$residuals)/sd(v_mod$residuals))

plot(h_mod$residuals ~ v_mod$residuals, asp=1)
plot(h_mod_sd ~ v_mod_sd)
plot(h_res ~ v_res, asp = 1, xlab = 'Vertical Residuals', ylab = 'Horzontal Residuals', main = 'Standardized Residuals')
plot(h_res ~ v_res, xlab = 'Vertical Residuals', ylab = 'Horzontal Residuals', main = 'Standardized Residuals',
     xlim = c(-3, 3), ylim = c(-3,3))
abline(v = -2)
abline(v = 2)
abline(h = -2)
abline(h = 2)

hres_95 <- ifelse(abs(h_res) < 2, h_res, 0)
vres_95 <- ifelse(abs(v_res) < 2, v_res, 0)
cor(h_res, v_res)
cor(hres_95, vres_95)

resids <- matrix(c(h_mod$residuals, v_mod$residuals), ncol = 2, byrow = F)

mshapiro.test(t(resids))
mvnorm.etest(resids, R = 100)
hoeffd(h_res, v_res)

summary(lm(h_res ~ v_res))

