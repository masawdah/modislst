
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fill the gaps in MODIS LST data

<!-- badges: start -->

<!-- badges: end -->

`modislst` is an R package for filling the gaps in MODIS LST data.
Different functions have been implemented to reconstruct Aqua/Terra LST
data using the nearest scene available in time, reconstruct the scene in
space using GAM with 3-dimensional spatial surface and additional
covariables like digital elevation model, or combine the two method
together to get better results.

## Installation

You can install the released version of modislst from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("modislst")
```

## Example

This is an example which shows you how to fill the gaps of daytime Aqua
scene (target) on 31 October 2016 using the available Terra scene in the
same day, then fill the remaining gaps completely in space.

### Setup

Call `modislst` and prepare the necessary data. In this example
Aqua/Terra scenes are projected and cropped already and matched with the
covariables (Digital Elevation Model and aspect).

``` r
library(modislst)
library(raster)
#> Loading required package: sp

## Add the target scene (aqua), and the reference scene (terra)
aqua <- raster("G:/Munster/Thesis/NewData/20161031/aquaDay.tif")
terra <- raster("G:/Munster/Thesis/NewData/20161031/terraDay.tif")

## Add covariables to improve the results 
dem <- raster("G:/Munster/Thesis/Data/DEM_resampled.tif")
aspect <- raster("G:/Munster/Thesis/Data/aspect.tif")

## Plot the Aqua/Terra scene before the filling process
plot(aqua,main="Aqua Scene Before Filling - Target Scene")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot(terra, main="Terra Scene - Refernce Scene")
```

<img src="man/figures/README-example-2.png" width="100%" />

### Fill the gaps in time

``` r
## Fill the gaps partially in time using moving window of size 47*47.
filledAqua <- fillingInTime(inputs = list(aqua,terra), m=47)

## Plot the filled Aqua scene 
plot(filledAqua, main="Filled Aqua Scene in time")
```

<img src="man/figures/README-reconstruct in time-1.png" width="100%" />

### Fill the gaps in space

``` r

## Fill the remaining gaps partially in space using.
filledAqua <- fillingInSpace(lstInput=filledAqua, covariable1=dem, covariable2=aspect)

## Plot the filled Aqua scene 
plot(filledAqua, main="Filled Aqua Scene")
```

<img src="man/figures/README-reconstruct in space-1.png" width="100%" />
