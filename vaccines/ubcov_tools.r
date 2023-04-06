###########################################################
### Author: USERNAME
### Date: DATE
### Project: ubCov
### Purpose: ubCov general utility functions
###########################################################

###################
### Setting up ####
###################
os <- .Platform$OS.type
if (os == "windows") {
  jpath <- "FILEPATH"
} else {
  jpath <- "FILEPATH"
}
package_lib <- paste0('FILEPATH')

####################################################################################################################################################
# 															   Table of Contents
####################################################################################################################################################

## Base
	## ubcov_path
	## ubcov_functions

####################################################################################################################################################
# 															         Base
####################################################################################################################################################


ubcov_path <- function(obj) {

	## Central Root 
	central_root <- paste0(unlist(strsplit(getwd(), "ubcov_central"))[1], "ubcov_central")

	if (obj == "central_root") {
		path_to_obj <- central_root
	} else {

		## Load paths file
		paths <- fread(paste0(central_root, "/paths.csv"))

		## Check that object is in path list
		if (length(paths[object==obj, object]) != 1) stop(paste0("Check that object: ", paste0(obj), " is in ~/paths.csv"))

		## Collapse paths together
		paths$path_full <- str_trim(apply(paths[,grep("path", names(paths), value=T), with=F], 1, paste, collapse=" "))
		
		## Find object
		path_to_obj <- as.character(paths[object==obj, .(path_full)])

		## Find parent object if there's a "," character
		while (grepl(",", path_to_obj)) {

			## Find path dependencies
			dependency <- unlist(strsplit(path_to_obj, ","))[1]

			## Run function to find path to dependency
			if (dependency != "central_root") {
				assign(dependency, ubcov_path(dependency))
			}

			## Replace with path
			path_to_obj <- paste0(get(dependency), "/", str_trim(unlist(strsplit(path_to_obj, ","))[2]))
		}
		}

	return(path_to_obj)
}

#####################################################################################################################################################

ubcov_functions <- function(list) {
	## Load functions
	for (func in list) {
		source(ubcov_path(func))
	}
}

#####################################################################################################################################################

path_loader <- function(file) {

	## Iterate through and set path locals into global environment
	## where full paths (root = FILEPATH) are set directly and recursive 
	## paths (ie blah_path = FILEPATH) are set through replacement

	## OS locals
    os <- .Platform$OS.type
    if (os == "windows") {
      jpath <- "FILEPATH"
    } else {
      jpath <- "FILEPATH"
    }

    ## Set paths
	paths <- fread(file)
	for (j in 2:ncol(paths)) {
	for (i in 1:nrow(paths)) {
		obj <- paths[i, 1, with=F] %>% unlist
		path <- paths[i, j, with=F] %>% unlist %>% gsub("FILEPATH", jpath, .)
		if (path != "" & !is.na(path)) assign(obj, path, envir=globalenv())
		if (grepl("[,]", path)) {
			replace <- gsub("[,].*$", "", path) %>% gsub(" ", "", .)
			end <- gsub(".*[,]", "", path) %>% gsub(" ", "", .)
			path <- paste0(get(replace), end)
			assign(obj, path, envir=globalenv())
		}
	}
	}

	## Done
	print("Paths loaded")

}
