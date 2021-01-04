

## code to prepare `DATASET` dataset goes here

path <- "./inst/"


## MODIS LST data
aqua <- raster("G:/Munster/Thesis/NewData/20161031/aquaDay.tif")
terra <- raster("G:/Munster/Thesis/NewData/20161031/terraDay.tif")

writeRaster(aqua, filename = file.path(paste0(path,"MODIS_data/"), "aqua.tif"))
writeRaster(terra, filename = file.path(paste0(path,"MODIS_data/"), "terra.tif"))



## Covariables
dem <- raster("G:/Munster/Thesis/Data/DEM_resampled.tif")
aspect <- raster("G:/Munster/Thesis/Data/aspect.tif")

raster::writeRaster(dem, filename = file.path(paste0(path,"covariables/"), "dem.tif"))
raster::writeRaster(aspect, filename = file.path(paste0(path,"covariables/"), "aspect.tif"))

