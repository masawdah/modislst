

#' @title Spatial reconstruction
#'
#' @description Fill the gaps in the MODIS LST scene using GAM with 3D spatial trend surface, and two covariables used  to improve the results. The covariables have to be raster objects and matching the input MODIS LST scene.
#'
#'
#' @param lstInput MODIS LST scene to be filled in space, using GAM with 3D spatial surface trend.
#' @param covariable1 Additional variable to improve the prediction of LST values. It should be a raster object and match the extent and origin of "lstInput" variable. Digital elevation model of Antarctica used as default value, because this package developed after analyzing MODIS LST data for Antarctica and found that DEM is very important to predict LST in Antarctica. The user can use DEM for his study area or even use another variable.
#'
#' @param covariable2 Additional variable to improve the prediction of LST values. It should be a raster object and match the extent and origin of "lstInput" variable. Aspect of Antarctica used as default value, because this package developed after analyzing MODIS LST data for Antarctica and found that aspect is very important to predict LST in Antarctica. The user can use aspect for his study area or even use another variable.
#'
#' @param nc Cluster size (Default = 4).
#' @param bs Represent the smooth function and quadratic penalty to control degree of freedom of smoothing (Default value is cubic spline basis "cr").
#' @param k  Number of knots used for smoothing the function and decide the degree of freedom (Default = 20).
#'
#'
#' @return Filled MODIS LST scene.
#' @export
#' @importFrom raster getData projection stack
#' @importFrom sp spTransform
#' @importFrom sf st_as_sf
#' @importFrom parallel detectCores makeCluster
#' @importFrom exactextractr exact_extract
#' @importFrom mgcv bam predict.bam
#' @importFrom stats complete.cases
#'
#'
#' @examples
#' \dontrun{
#' fillingInSpace(aqua, covariable1=dem, covariable1=aspect)
#' }
#'
#'
#'

fillingInSpace <- function(lstInput, covariable1=dem, covariable2=aspect, nc=4, bs="cr", k=20){
  ##binding the "cl" variable locally to the function
  cl <- NULL
  ## Check the number of available cores for parallel computing, and prepare the clusters for parallel computing
  funs <- function(nc=4){

    if (detectCores() > nc) { ## no point otherwise
      message(paste0("Your machine has ",detectCores(), " cores, ","you can use all of them instead of ", nc ," cores to decrease the processing time"))
      cl <- makeCluster(nc)

    } else if(detectCores() == nc) {
      cl <- makeCluster(nc)
    } else {
      stop(paste0("Your machine has ",detectCores(), " cores, ","please use the same number of cores or less"))
    }

  }


  ## Match all variables in one stack
  names(lstInput) <- c("lst")
  names(covariable1) <- c("ele")
  names(covariable2) <- c("aspect")
  s <- stack(lstInput,covariable1,covariable2)

  ## convert the polygon to simple feature geometry (to make the extraction values faster)
  studyarea <- getData(country = "ATA", level = 0, path = system.file("/inst/studyarea",package="modislst"))
  studyarea <- spTransform(studyarea,projection(lstInput))
  studyarea <- st_as_sf(studyarea)

  ## extract LST values
  dataset <- exact_extract(s, studyarea,include_cell=TRUE,include_xy = TRUE)
  dataset <- dataset[[1]]

  ## remove pixels not covered by digital elevation model
  dataset <- dataset[complete.cases(dataset[ 2:5]),]

  ## extract the missing pixels
  missing.pixels <- dataset[is.na(dataset[ , 1]),]

  ## Full dataset (All pixels with valid LST values)
  dataset <- dataset[complete.cases(dataset[ ,1]),]


  ## Build GAM with 3D spatial trend
  model<- bam(lst ~ s(ele,bs=bs,k=k)+aspect+te(x,y,ele,bs=c("tp","tp","cr")),method ="REML",data=dataset,cluster=cl)

  ## Predict values of missing pixels
  predictions <- predict.bam(model,missing.pixels,cluster=cl)
  lstInput[missing.pixels$cell] <- predictions

  ## Return filled MODIS LST scene
  return(lstInput)



}
