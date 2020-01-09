# Joshua North
# 1/23/2017


################################################################
### Move Functions
################################################################

# obs1 is the initial/first image/timestep you wish to use
# obs2 is the second timestep you wish to use
# nx is the dimension in the x direction
# ny is the dimension in the y direction
# size is the fraction of pixel you want to move, for example size = 1 is a whole pixel and size = .1 is 1/10 of a pixel.
# bigD is the array of [nx,ny,number of images] where your data is stored, must be a 3-D matrix of the format: "matrix"    "array"     "structure" "vector" 
# tolerance is the amount of error you wish to allow, setting this to zero ensures that it will run until iteration has hit its mark
# number is the fraction you want your accuracy to be to, for example, number = 50 means the smallest movement could be 1/50 of a pixel

################################################################
### Move Functions
################################################################

# Shifts matrix up and image left one unit
m.left <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mleft <- matrix(c(0),nr = nx, nc = ny) # initialize matrix
	mleft[1:(nx-1),] <- ((obs1[-1,]*(size)) + (obs1[-nx,]*(1-size))) # take a portion of left image plus portion of image and add
	mleft[nx,] <- obs1[nx,] # interpolate by adding in far right column
	error <- ERROR(mleft,obs2) # calculate the error
	list(mleft = mleft, error = error) # output matrix and error
}

# Shifts matrix down and image right one unit
m.right <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mright <- matrix(c(0),nr = nx, nc = ny)
	mright[2:nx,] <- ((obs1[-nx,]*(size)) + (obs1[-1,]*(1-size)))
	mright[1,] <- obs1[1,]
	error <- ERROR(mright,obs2)
	list(mright = mright,  error = error)
}

# Shifts matrix right and image up one unit
m.plus <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mplus <- matrix(c(0),nr = nx, nc = ny)
	mplus[,2:ny] <- ((obs1[,-ny]*(size)) + (obs1[,-1]*(1-size)))
	mplus[,1] <- obs1[,1]
	error <- ERROR(mplus,obs2)
	list(mplus = mplus, error = error)
}

# Shifts matrix left and image down one unit
m.minus <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mminus <- matrix(c(0),nr = nx, nc = ny)
	mminus[,1:(ny-1)] <- ((obs1[,-1]*(size)) + (obs1[,-ny]*(1-size)))
	mminus[,-ny] <- mminus[,-ny]
	error <- ERROR(mminus,obs2)
	list(mminus = mminus, error = error)
}

# Shifts matrix up/right and image left/up
m.ul <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mul <- matrix(c(0),nr = nx, nc = ny)
	mul[1:(nx-1),2:ny] <- ((obs1[-1,-ny]*(size)) + (obs1[-nx,-1]*(1-size)))
	mul[nx,] <- obs1[nx,]
	mul[,1] <- obs1[,1]
	error <- ERROR(mul,obs2)
	list(mul = mul, error = error)
}

# Shifts matrix up/left and image left/down
m.dl <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mdl <- matrix(c(0),nr = nx, nc = ny)
	mdl[1:(nx-1),1:(ny-1)] <- ((obs1[-1,-1]*(size)) + (obs1[-nx,-ny]*(1-size)))
	mdl[nx,] <- obs1[nx,]
	mdl[,ny] <- obs1[,ny]
	error <- ERROR(mdl,obs2)
	list(mdl = mdl, error = error)
}

# Shifts matrix down/right and image right/up
m.ur <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mur <- matrix(c(0),nr = nx, nc = ny)
	mur[2:nx,2:ny] <- ((obs1[-nx,-ny]*(size)) + (obs1[-1,-1]*(1-size)))
	mur[1,] <- obs1[1,]
	mur[,1] <- obs1[,1]
	error <- ERROR(mur,obs2)
	list(mur = mur, error = error)
}

# Shifts matrix down/left and image right/down
m.dr <- function(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = 1){
	mdr <- matrix(c(0),nr = nx, nc = ny)
	mdr[2:nx,1:(ny-1)] <- ((obs1[-nx,-1]*(size)) + (obs1[-1,-ny]*(1-size)))
	mdr[1,] <- obs1[1,]
	mdr[,ny] <- obs1[,ny]
	error <- ERROR(mdr,obs2)
	list(mdr = mdr, error = error)
}


################################################################
### Move Whole Image Function
################################################################

move.all <- function(obs1, obs2, tolerance = 0, iteration = 10, number = 10){

	counter = 0			# initialize counter to zero
	tmp <- dim(obs1)		# temp variable for dimensions
	nx <- tmp[1]		# x dimension
	ny <- tmp[2]		# y dimension
	error.orig <- error <- ERROR(obs1, obs2)		# original error between obs1 and obs2
	move <- matrix(c(0), nc = 2, nr = 1)		# initialize the movement matrix
	colnames(move) <- c("Horizontal", "Vertical")	# add column names


	while((error >= tolerance) && (counter < iteration)){		# while error is less than tolerance and counter is less than iteration


		# Note that all of the movements are with respect to the image not the matrix.
		# The image is indexed with [1,1] being the bottom left corner while the matrix
		# is indexed with [1,1] being the top left corner.

		error <- rep(0, number*8)				# allocate array for error to be stored in
		matrix <- array(NA, dim=c(nx, ny, number*8))	# allocate 3D matrix to store moved matrices in, there is enough space for each direction at each fractional value of the move
		tmp <- matrix(c(0), nc = 2, nr = number*8)	# allocate array for the amount each matrix move value goes

		# Move Left
		for(i in 1:number){	# for each of the fractional movements, i.e if number is 10 there will be a move of 1/10, 2/10,...,10/10
			mleft <- m.left(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)		# call m.left
			error[i] <- mleft$error			# assign error to error array
			matrix[,,i] <- mleft$mleft		# assign the moved matrix to the 3D matrix store
			tmp[i,1] <- tmp[i,1] - i/number	# amount moved, minus because it is left, plus is right
		}
		
		for(i in 1:number){
			mright <- m.right(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)
			error[(i+number)] <- mright$error
			matrix[,,(i+number)] <- mright$mright
			tmp[(i+number),1] <- tmp[(i+number),1] + i/number
		}

		for(i in 1:number){
			mminus <- m.minus(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)
			error[(i+2*number)] <- mminus$error
			matrix[,,(i+2*number)] <- mminus$mminus
			tmp[(i+2*number),2] <- tmp[(i+2*number),2] - i/number
		}
		
		for(i in 1:number){
			mplus <- m.plus(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)
			error[(i+3*number)] <- mplus$error
			matrix[,,(i+3*number)] <- mplus$mplus
			tmp[(i+3*number),2] <- tmp[(i+3*number),2] + i/number
		}

		for(i in 1:number){
			mul <- m.ul(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)
			error[(i+4*number)] <- mul$error
			matrix[,,(i+4*number)] <- mul$mul
			tmp[(i+4*number),1] <- tmp[(i+4*number),1] - i/number
			tmp[(i+4*number),2] <- tmp[(i+4*number),2] + i/number
		}
		
		for(i in 1:number){
			mdl <- m.dl(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)
			error[(i+5*number)] <- mdl$error
			matrix[,,(i+5*number)] <- mdl$mdl
			tmp[(i+5*number),1] <- tmp[(i+5*number),1] - i/number
			tmp[(i+5*number),2] <- tmp[(i+5*number),2] - i/number
		}

		for(i in 1:number){
			mur <- m.ur(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)
			error[(i+6*number)] <- mur$error
			matrix[,,(i+6*number)] <- mur$mur
			tmp[(i+6*number),1] <- tmp[(i+6*number),1] + i/number
			tmp[(i+6*number),2] <- tmp[(i+6*number),2] + i/number
		}
		
		for(i in 1:number){
			mdr <- m.dr(obs1 = obs1, obs2 = obs2, nx = nx, ny = ny, size = i/number)
			error[(i+7*number)] <- mdr$error
			matrix[,,(i+7*number)] <- mdr$mdr
			tmp[(i+7*number),1] <- tmp[(i+7*number),1] + i/number
			tmp[(i+7*number),2] <- tmp[(i+7*number),2] - i/number
		}

		temp <- order(error)		# order based on least amount of error on top, allows you select values easier now
		obs1.new <- matrix[,,temp[1]] # select the moved matrix that correlates to the smallest error
		move.new <- tmp[temp[1],]	# select the amount moved that correlates to the smallest error
		error.new <- error[temp[1]]	# select the amount for smallest error

		if( error.orig >= error.new){	# if the new error is less than original error do
			obs1 <- obs1.new		# assign obs1 to get the new obs for the continuation of the loop
			error <- error.new	# assign error to get new error
			move <- move + move.new	# move gets move plus new move value
		}else{				# if original error is less than new error break
			break
		} # end if

		counter = counter +1		# incriment counter
	} # end while

	error <- ERROR(obs1,obs2)		# calculate new error between moved obs1 and original obs2
	list(obs1 = obs1, move = move, error = error)	# output moved obs1, amount moved and error


} # end move.all


################################################################
### Show where points move
################################################################

# Shows you how far the image moves for as many time steps as you want
points.move.partial <- function(bigD = bigD, tolerance = 0, iteration = 10, number = 10, timestep = timestep){

	store <- matrix(c(0), nr = (timestep-1), nc = 3)
	colnames(store) <- c("Horizontal", "Vertical", "Error")
	for (i in 1: (timestep-1)){
		obs1 <- bigD[ , ,i]
		obs2 <- bigD[ , ,(i+1)]
		tmp <- move.all(obs1, obs2, tolerance = tolerance, iteration = iteration, number = number )
		store[i,1] <- tmp$move[1]
		store[i,2] <- tmp$move[2]
		store[i,3] <- tmp$error
	}
	store<- cbind(store,"order"=1:nrow(store))
	return(store)

}


