

#' @title Spatial reconstruction
#'
#' @description Fill the MODIS LST scene gaps using GAM with 3D spatial trend surface and two covariables to improve the results. The covariables have to be raster objects and matching the input MODIS LST scene.
#'
#'
#' @param lstInput MODIS LST scene to be filled in space.
#' @param studyarea Name of the study area ( 'GADM','countries'). The default study area is Antarctica "ATA".
#' @param covariable1 Digital elevation model as an additional variable to build the 3D spatial trend surface. It should be a raster object and match the extent and origin of the "lstInput" variable.
#'
#' @param covariable2 Any additional variable to improve the prediction of LST values. It should be a raster object and match the extent and origin of the "lstInput" variable.
#'
#' @param nc Cluster size (Default = 4).
#' @param bs Represent the smooth function and quadratic penalty to control the degree of freedom of smoothing (Default value is cubic spline basis "cr").
#' @param k  The Number of knots used for smoothing the function and decide the degree of freedom (Default = 20).
#'
#'
#' @return Filled MODIS LST scene.
#' @export
#' @importFrom raster getData projection stack raster
#' @importFrom sp spTransform
#' @importFrom sf st_as_sf
#' @importFrom parallel detectCores makeCluster
#' @importFrom exactextractr exact_extract
#' @importFrom mgcv bam  predict.bam
#' @importFrom stats complete.cases
#'
#'
#' @examples
#'\dontrun{
#' aqua <- raster(system.file("MODIS_data","aqua.tif", package="modislst"))
#' dem <- raster(system.file("covariables","dem.tif", package="modislst"))
#' aspect <- raster(system.file("covariables","aspect.tif", package="modislst"))
#' fillingInSpace(aqua, covariable1=dem, covariable1=aspect)
#'}
#'
#'
#'

fillingInSpace <- function(lstInput, studyarea = "ATA", covariable1, covariable2, nc=4, bs="cr", k=20){
  ##binding the "cl" variable locally to the function
  cl <- NULL

  ## Check the number of available cores for parallel computing, and prepare the clusters for parallel computing
  if (detectCores() > nc) {
    message(paste0("Your machine has ",detectCores(), " cores, ","you can use all of them instead of ", nc ," cores to decrease the processing time"))
    cl <- makeCluster(nc)
    } else if(detectCores() == nc) {
      cl <- makeCluster(nc)
      } else {
        stop(paste0("Your machine has ",detectCores(), " cores, ","please use the same number of cores or less"))
        }

  ## Match all variables in one stack
  names(lstInput) <- c("lst")
  names(covariable1) <- c("ele")
  names(covariable2) <- c("cov2")
  s <- stack(lstInput,covariable1,covariable2)

  ## convert the polygon to simple feature geometry (to make the extraction values faster)
  studyarea <- getData(country = studyarea, level = 0, path = "inst/studyarea/")
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
  model<- bam(lst ~ s(ele,bs=bs,k=k)+s(cov2,bs=bs,k=k)+te(x,y,ele,bs=c("tp","tp","cr")),method ="REML",data=dataset,cluster=cl)

  ## Predict values of missing pixels
  predictions <- predict.bam(model,missing.pixels,cluster=cl)
  lstInput[missing.pixels$cell] <- predictions

  ## Return filled MODIS LST scene
  return(lstInput)



}
