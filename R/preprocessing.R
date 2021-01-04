
#' @title Pre-processing function
#'
#'
#' @description        Crop MODIS LST scenes based on the extent of the study area,after that reproject cropped MODIS LST scenes to the target projected coordinate system.
#'
#' @param targetLST    The input MODIS LST scene to be filled, downloaded from MODIS package.
#' @param referenceLST The input MODIS LST scene to fill the gaps from it, downloaded from MODIS package.
#' @param studyarea    Name of study area ( 'GADM','countries'), the default study area is Antarctica "ATA".
#' @param crs          Target projected coordinate system, the default is Antarctica polar stereographic projected coordinate system "+init=epsg:3031".
#'
#' @return List of cropped and projected MODIS LST scenes, the first element in the list is the target scene and the second one is the reference scene.
#' @export
#' @importFrom raster getData projection crop projectRaster
#' @importFrom sp spTransform
#'
#' @examples
#' \dontrun{
#' lstProcessing(aqua,terra)
#' }
#'
#'
lstProcessing <- function(targetLST, referenceLST, studyarea="ATA", crs="+init=epsg:3031"){
  studyarea <- getData(country = studyarea, level = 0,path = "inst/studyarea/")
  studyarea <- spTransform(studyarea,projection(targetLST))


  ## crop LST scene using study area extent
  targetLST <- crop(targetLST,studyarea)
  referenceLST <- crop(referenceLST,studyarea)

  ## Re-project LST scene
  targetLST <- projectRaster(targetLST,res = 1000, crs=crs(crs))
  referenceLST <- projectRaster(referenceLST,res = 1000, crs=crs(crs))


  lst <- list(targetLST, referenceLST)
  return(lst)
}
