# Joshua North
# 1/23/2017

################################################################
### Load netCDF Data Function
################################################################

# filepath is the file path of where the data is stored.  
# Must enter as dummy <- load.netCDF.data(filepath = 'filepath')
# Must be netCDF data

################################################################
### Load netCDF Data Function
################################################################

load.netCDF.data <- function(filepath = filepath){

	###### Load in data, data from Eric is in 5 min incriments 
	files= list.files(filepath,pattern='*.nc',full.names=TRUE)

	###### Empty matrix to fill in
	store <- array(NA, dim=c(400, 400, 72))

	###### fill in matrix
	for(i in seq_along(files)) {
    	nc <- open.nc(files[i], write = FALSE)
	    store[,,i] <- var.get.nc(nc, "LtgPotentialclu")
    	store[,,i][is.na(store[,,i])] <- 0
    	close.nc(nc)
		rm(nc)
	}
	
	return(store)

}