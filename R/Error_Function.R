# Joshua North
# 1/23/2017

################################################################
### Error Function
################################################################

# obs1 is the initial/first image/timestep you wish to use
# obs2 is the second timestep you wish to use

################################################################
### Error Function
################################################################

ERROR <- function(obs1, obs2){
	error <- sqrt(mean((obs2 - obs1)^2))
	return(error)
}