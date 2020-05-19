# Joshua North
# 1/23/2017

################################################################
### Prediction Functions
################################################################

# obs1 is the initial/first image/timestep you wish to use
# obs2 is the second timestep you wish to use
# nx is the dimension in the x direction
# ny is the dimension in the y direction
# size is the fraction of pixel you want to move, for example size = 1 is a whole pixel and size = .1 is 1/10 of a pixel.
# bigD is the array of [nx,ny,number of images] where your data is stored, must be a 3-D matrix of the format: "matrix"    "array"     "structure" "vector" 
# tolerance is the amount of error you wish to allow, setting this to zero ensures that it will run until iteration has hit its mark
# iteration is the how far you think the two images are appart pixel-wise, for instance, iteration = 10 means you think two images will be at most 10 pixels appart/shifted
# number is the fraction you want your accuracy to be to, for example, number = 50 means the smallest movement could be 1/50 of a pixel
# timestep is the timepoint where your last image is, so if you choose timestep = 50 it means you will go up to bigD[,,50] then predict your images from there
# step is how many steps in the future you want to go, if step = 10 and timestep = 50 you will predict to the 60th timestep, or 10 x time increment between images


# Prediction function based off of a regression
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
	
	# mod <- lm(cbind(Horizontal, Vertical) ~ order, data = run)

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
	  
	  # preds <- predict(mod, newdata = data.frame(order = (timestep+k)))
	  # horizontal = preds[,1]
	  # vertical = preds[,2]
	  
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



# Prediction function based off of the amount moved from the last image
predict.sustain <- function(bigD = bigD, tolerance = 0, iteration = 10, number = 10, 
                            timestep = timestep, step = 10){

	run <- points.move.partial(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = timestep)
	a <- run[,1] # horizontal movement
	b <- run[,2] # vertical movement
	c <- run[,3] # error
	d <- run[,4] # order

	hor <- lm(a~d)
	ver <- lm(b~d)
	
	curr <- bigD[,,timestep]
	tmp <- dim(curr)
	nx <- tmp[1]
	ny <- tmp[2]

	store <- array(NA, dim=c(400,400,step))
	
	for ( k in 1:step){

		if(k == 1){
			curr <- bigD[,,timestep]
		}else{
			curr <- store[,,k-1]
			}
		
		horizontal = run[(timestep-1), 1]
		vertical = run[(timestep-1), 2]
		r_hor <- floor(horizontal)
		r_ver <- floor(vertical)
		hor_remainder <- horizontal - r_hor
		ver_remainder <- vertical - r_ver
		
		nex <- matrix(c(0), nr = nx, nc = ny)

		if(horizontal > 0){
			if(vertical > 0){
					
				for(i in 1:(nx - r_hor) ){ 			# move integer value
					for(j in 1:(ny - r_ver) ){
						
						nex[i,j] <- curr[(i+r_hor),(j+r_ver)]

					} # end for j
				} # end for i 
				
				nex[(1:r_hor),(1:r_ver)] <- curr[(1:r_hor),(1:r_ver)] 			# interpolate
				tmp <- m.right(obs1 = nex, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
				tmp_2 <- m.plus(obs1 = tmp$mright, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
				store[,,k] <- tmp_2$mplus 		# put into store matrix
				
			}else{ # if vertical < 0, end vertical if
					
				for(i in 1:(nx - r_hor) ){
					for(j in (1 + r_ver):ny ){
						
						nex[i,j] <- curr[(i+r_hor),(j-r_ver)]
						
					} # end for j
				} # end for i	
				
				nex[(1:r_hor),((ny - r_ver) : ny)] <- curr[(1:r_hor),((ny - r_ver) : ny)]			# interpolate
				tmp <- m.right(obs1 = nex, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
				tmp_2 <- m.minus(obs1 = tmp$mright, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
				store[,,k] <- tmp_2$mminus 		# put into store matrix
				
			} # end else 
			
		}else{ # horizontal < 0 end horizontal if
			if(vertical > 0){
					
				for(i in (1 - r_hor):nx){
					for(j in 1:(ny - r_ver) ){
						
						nex[i,j] <- curr[(i-r_hor),(j+r_ver)]
						
						} # end j for
					} # end i for
									
				nex[((nx - r_hor) : nx),(1 : r_ver)] <- curr[((nx - r_hor) : nx),(1 : r_ver)]  			# interpolate
				tmp <- m.left(obs1 = nex, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
				tmp_2 <- m.plus(obs1 = tmp$left, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
				store[,,k] <- tmp_2$mplus 		# put into store matrix
					
			}else{ #vertical < 0
					
				for(i in (1-r_hor):nx){
					for(j in (1-r_ver):ny){
						
						nex[i,j] <- curr[(i-r_hor),(j-r_ver)]
				
					} # end j for
				} # end i for
				
				nex[((nx - r_hor) : nx),((ny - r_ver) : ny)] <- curr[((nx - r_hor) : nx),((ny - r_ver) : ny)]  			# interpolate
				tmp <- m.left(obs1 = nex, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
				tmp_2 <- m.minus(obs1 = tmp$left, obs2 = bigD[,,timestep],nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
				store[,,k] <- tmp_2$mminus 		# put into store matrix

			} # end vertical < 0 else
		} # end horizontal < 0 else
	}
	
	return(store)
	
} # end function




# error_analysis_function -------------------------------------------------

# Prediction function based off of a regression
predict.regression.error <- function(obs, tolerance = 0, iteration = 10, number = 10, timestep, step = 10, mu, sigma){
    
  
    curr <- obs
    
    tmp <- dim(curr)
    nx <- tmp[1]
    ny <- tmp[2]
    
    store <- array(NA, dim=c(400,400,step))
    
    vals = mvrnorm(1, mu, sigma)
    
    for ( k in 1:step){
      
      horizontal = vals[1]
      vertical = vals[2]
      if(horizontal > 400){
        horizontal <- 399
      }else if(horizontal < -400){
        horizontal <- -399
      }
      if(vertical > 400){
        vertical <- 399
      }else if(vertical < -400){
        vertical <- -399
      }
      r_hor <- floor(horizontal)
      r_ver <- floor(vertical)
      hor_remainder <- horizontal - r_hor
      ver_remainder <- vertical - r_ver
      
      nex <- matrix(c(0), nr = nx, nc = ny)
      
      if(horizontal > 0){		# to move image right
        if(vertical > 0){		# to move image up
          
          for(i in 1 : (nx - r_hor) ){			# choose 1 : nx-r_hor (whole number) of matrix
            for(j in (1 + r_ver) : ny ){		# choose 1+r_ver (whole number) of matrix
              
              nex[j,i] <- curr[(j-r_ver),(i+r_hor)]		# shift matrix down and right (moves image up and right)
              
            } # end for j
          } # end for i	
          
          nex[(1:r_hor),((ny - r_ver) : ny)] <- curr[(1:r_hor),((ny - r_ver) : ny)]			# interpolate
          tmp <- m.right(obs1 = nex, obs2 = curr,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
          tmp_2 <- m.plus(obs1 = tmp$mright, obs2 = curr,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
          store[,,k] <- tmp_2$mplus 		# put into store matrix
          
        }else{ # if vertical < 0, end vertical if
          
          for(i in 1 : (nx - r_hor) ){ 			# move integer value
            for(j in 1 : (ny + r_ver) ){
              
              nex[j,i] <- curr[(j-r_ver),(i+r_hor)]
              
            } # end for j
          } # end for i 
          
          nex[(1:r_hor),(1:abs(r_ver))] <- curr[(1:r_hor),(1:abs(r_ver))] 			# interpolate
          tmp <- m.right(obs1 = nex, obs2 = curr,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
          tmp_2 <- m.minus(obs1 = tmp$mright, obs2 = curr,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
          store[,,k] <- tmp_2$mminus 		# put into store matrix
          
        } # end else 
        
      }else{ # horizontal < 0 end horizontal if
        if(vertical > 0){
          
          for(i in (1 - r_hor) : nx){
            for(j in (1 + r_ver) : ny){
              
              nex[j,i] <- curr[(j-r_ver),(i+r_hor)]
              
            } # end j for
          } # end i for
          
          nex[((nx - abs(r_hor)) : nx),((ny - r_ver) : ny)] <- curr[((nx - abs(r_hor)) : nx),((ny - r_ver) : ny)]  			# interpolate
          tmp <- m.left(obs1 = nex, obs2 = curr ,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
          tmp_2 <- m.plus(obs1 = tmp$mleft, obs2 = curr , nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
          store[,,k] <- tmp_2$mplus 		# put into store matrix
          
        }else{ #vertical < 0
          
          for(i in (1 + abs(r_hor)) : nx){
            for(j in 1 : (ny - abs(r_ver)) ){
              
              nex[j,i] <- curr[(j+abs(r_ver)),(i-abs(r_hor))]
              
            } # end j for
          } # end i for
          
          nex[((nx - abs(r_hor)) : nx),(1 : abs(r_ver))] <- curr[((nx - abs(r_hor)) : nx),(1 : abs(r_ver))]  			# interpolate
          tmp <- m.left(obs1 = nex, obs2 = curr ,nx = nx, ny = ny, size = hor_remainder)			 # move fractional value
          tmp_2 <- m.minus(obs1 = tmp$mleft, obs2 = curr ,nx = nx, ny = ny, size = ver_remainder) 			# move fractional value
          store[,,k] <- tmp_2$mminus 		# put into store matrix
          
        } # end vertical < 0 else
      } # end horizontal < 0 else
      
      curr <- store[,,k]
    }
    
    list(mat = store, mu1 = horizontal, mu2 = vertical)
    
  } # end function
  
