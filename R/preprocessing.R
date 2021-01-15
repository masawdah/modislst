
#' @title Pre-processing function
#'
#'
#' @description        Crop MODIS LST scenes based on the extent of the study area,after that reproject cropped MODIS LST scenes to the target projected coordinate system.
#'
#' @param LSTScene     MODIS LST scene to be cropped and projected.
#' @param studyarea    Name of study area ( 'GADM','countries'), the default study area is Antarctica "ATA".
#' @param crs          Target projected coordinate system, the default is Antarctica polar stereographic projected coordinate system "+init=epsg:3031".
#'
#' @return  Projected and Cropped MODIS LST scenes.
#' @export
#' @importFrom raster getData projection crop projectRaster
#' @importFrom sp spTransform
#'
#' @examples
#'
#' lstProcessing(aqua,terra)
#'
#'
#'
lstProcessing <- function(LSTScene, studyarea="ATA", crs="+init=epsg:3031"){
  studyarea <- getData(country = studyarea, level = 0,path = "inst/studyarea/")
  studyarea <- spTransform(studyarea,projection(LSTScene))


  ## crop LST scene using study area extent
  LSTScene <- crop(LSTScene,studyarea)


  ## Re-project LST scene
  LSTScene <- projectRaster(LSTScene,res = 1000, crs=crs(crs))


  return(LSTScene)
}
