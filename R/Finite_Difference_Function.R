# Joshua North
# 1/23/2017

################################################################
### Finite Difference Function
################################################################

# uInit is the last timestep of your data set, or the point you wish to predict from, used as initial image
# deltaT/deltaX needs to be less than 1/2 in 1D and 1/4 in 2D
# nframes is how many times the finite difference smooth is applied, see for loop
# Q is a smoothing matrix

# Finite Difference Function ----------------------------------------------

finiteD <- function(r = 0.2375, deltaT = .000001, deltaX = .0025, uInit = uInit, nframes = 20, Q){
#finiteD <- function(r = 0.0004, uInit = uInit, nframes = 20, Q){
  
  if(r){
    r = r
  }else{
    r = deltaT/deltaX^2 # must be less than 1/4 to be stable
  }
	
	if(r > 1/4) stop("must choose deltaT and deltaX such that deltaT/deltaX^2 < 1/4")

	temp = dim(uInit)
	nx = temp[1]
	ny = temp[2]
	uFut = matrix(c(0),nr = nx, nc = ny)

	for(k in 1 : nframes){
	
		for(i in 1 : ny){
		
			for(j in 1 : nx){
			
				if(i == 1 & j == 1){
				
					uFut[i,j] = uInit[i,j] + r*(2*uInit[i,j] + uInit[i+2,j] + uInit[i,j+2] - 2*uInit[i+1,j] - 2*uInit[i,j+1]) + Q[i,j]
				
				}else if(i == 1 & j == nx){
				
					uFut[i,j] = uInit[i,j] + r*(2*uInit[i,j] - 2*uInit[i+1,j] - 2*uInit[i,j-1] + uInit[i+2,j] + uInit[i,j-2]) + Q[i,j]
				
				}else if(i == ny & j == 1){
				
					uFut[i,j] = uInit[i,j] + r*(2*uInit[i,j] - 2*uInit[i,j+1] - 2*uInit[i-1,j] + uInit[i,j+2] + uInit[i-2,j]) + Q[i,j]
				
				}else if(i == ny & j == nx){
				
					uFut[i,j] = uInit[i,j] + r*(2*uInit[i,j] + uInit[i-2,j] + uInit[i,j-2] - 2*uInit[i-1,j] - 2*uInit[i,j-2]) + Q[i,j]
				
				}else if (j == 1){
				
					uFut[i,j] = uInit[i,j] + r*(uInit[i-1,j] - uInit[i,j] + uInit[i+1,j] - 2*uInit[i,j+1] + uInit[i,j+2]) + Q[i,j]
				
				}else if(j == nx){
				
					uFut[i,j] = uInit[i,j] + r*(uInit[i-1,j] - uInit[i,j] + uInit[i+1,j] - 2*uInit[i,j-1] + uInit[i,j-2]) + Q[i,j]
				
				}else if(i == 1){
				
					uFut[i,j] = uInit[i,j] + r*(uInit[i,j-1] - uInit[i,j] + uInit[i,j+1] - 2*uInit[i+1,j] + uInit[i+2,j]) + Q[i,j]
				
				}else if(i == ny){
				
					uFut[i,j] = uInit[i,j] + r*(uInit[i,j-1]  - uInit[i,j] + uInit[i,j+1] - 2*uInit[i-1,j] + uInit[i-2,j]) + Q[i,j]
				
				}else{
				
					uFut[i,j] = uInit[i,j] + r*(uInit[i+1,j]  + uInit[i-1,j]+ uInit[i,j+1] + uInit[i,j-1] - 4*uInit[i,j]) + Q[i,j]
				
				} # end if
			
			} # end for j
		
		} # end for i
		
		uInit <- uFut
		
	} # end for k

	return(uFut)

}

