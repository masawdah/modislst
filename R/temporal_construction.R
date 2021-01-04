

#' @title Temporal reconstruction
#'
#' @description Fill the gaps in the target MODIS LST scene from the nearest MODIS LST scene (reference scene).
#'
#' @param inputs List of length two with following inputs: the first element is the target scene, the second element is the reference scene.
#' @param m The size of moving window. Default value is 47, using bigger size would fill more pixels but increase the computing work and reduce the accuracy of filled pixels.
#'
#'
#' @return MODIS LST scene filling partially.
#' @export
#' @importFrom raster focal overlay raster
#' @importFrom stats quantile
#'
#' @examples
#'\dontrun{
#' aqua <- raster(system.file("MODIS_data","aqua.tif", package="modislst"))
#' terra <- raster(system.file("MODIS_data","terra.tif", package="modislst"))
#' fillingLST(inputs=list(aqua,terra), m=47)
#' }
#'
fillingInTime <- function(inputs,m=47){

  targetLST = inputs[[1]]
  referenceLST = inputs[[2]]

  ## calculate the differences between valid pixels in the two scenes
  lst.diff <- targetLST - referenceLST

  ## Exclude the extreme differences between the two scenes using quantile outlier detector
  qa <- quantile(lst.diff,na.rm=TRUE)
  IQR_qa <- qa[4] - qa[2]
  a <- IQR_qa*1.5
  b <- qa[4] + a
  c <- qa[2] - a
  lst.diff[lst.diff[] < c] <- NA
  lst.diff[lst.diff[] > b] <- NA

  local.mean.diff <- focal(lst.diff,w=matrix(1,m,m),na.rm=TRUE,pad=TRUE,fun=mean)

  fill.raster <- overlay(targetLST, referenceLST, local.mean.diff, fun = function(x, y, z) {
    x[is.na(x)] <- y[is.na(x[])] + z[is.na(x[])]
    return(x)
  })

  fill.raster <- (fill.raster*0.02) - 273.15



  return(fill.raster)

}
