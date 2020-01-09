# Joshua North
# 1/23/2017


####################################################################
### Nowcasting Function
####################################################################

## bigD is the array of images being looked at
## tolerance is there allowable error between two timesteps, the lower the tolerance the better the prediction
## iteration is the how far you think the two images are appart pixel-wise, for instance, iteration = 10 means you think two images will be at most 10 pixels appart/shifted
## number is the fraction you want your accuracy to be to, for example, number = 50 means the smallest movement could be 1/50 of a pixel
## timestep is the timepoint where your last image is, so if you choose timestep = 50 it means you will go up to bigD[,,50] then predict your images from there
## num is how many steps in the future you want to go, if step = 10 and timestep = 50 you will predict to the 60th timestep, or 10 x time increment between images
## deltaT/deltaX needs to be less than 1/2 in 1D and 1/4 in 2D
## uInit is the last timestep of your data set, or the point you wish to predict from

####################################################################
### Nowcasting Function
####################################################################

nowcast <- function(bigD = bigD, tolerance = 0, iteration = 10, number = 10,
                    timestep = timestep, num = 10, r = 0.2375, deltaT = .000001,
                    deltaX = .0025, nframes = 20, Q){

	store <- array(NA, dim=c(400,400,num))
	prev = 0

	dt = deltaT
	dx = deltaX
	q = Q
	nf = nframes


	for(i in 1 : num){

	  smpl1 <- predict.regression(bigD = bigD, tolerance = 0, iteration = 10, number = 10,
	                              timestep = timestep, prev = prev, step = 1, run = run)

		temp <- smpl1[,,1]
		smpl2 <- finiteD(deltaT = dt, deltaX = dx, uInit = temp, nframes = nf, Q = q)

		store[,,i] <- smpl2
		prev <- smpl2

	} # end for

	return(store = store)

} # end nowcast


# error bounds nowcast function
nowcast.error <- function(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = timestep, num = 10, deltaT = .000001, deltaX = .0025, nframes = 20, hor, ver){
  
  store <- array(NA, dim=c(400,400,num))
  prev = 0
  
  for(i in 1 : num){
    
    smpl1 <- predict.regression.error(bigD, tolerance, iteration, number, timestep, prev, step = 1, hor, ver)
    temp <- smpl1[,,1]
    smpl2 <- finiteD(deltaT, deltaX, uInit = temp, nframes)
    
    store[,,i] <- smpl2
    prev <- smpl2
    
  } # end for
  
  return(store = store)
  
} # end nowcast

