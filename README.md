
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Fill the gaps in MODIS LST data

<!-- badges: start -->

[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/masawdah/modislst?branch=master&svg=true)](https://ci.appveyor.com/project/masawdah/modislst)
[![Travis build
status](https://travis-ci.com/masawdah/modislst.svg?branch=master)](https://travis-ci.com/masawdah/modislst)
<!-- badges: end -->

`modislst` is an R package for filling the gaps in MODIS LST data.
Different functions have been implemented to reconstruct Aqua/Terra LST
data using the nearest scene available in time, reconstruct the scene in
space using GAM with 3-dimensional spatial surface and additional
covariables like digital elevation model and aspect, or combine the two
method together to get better results.

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

Call `modislst` and prepare the necessary data. This example will show
you how to fill Aqua LST scene for Antarctica.

In this example Aqua/Terra LST overpasses on 31 October 2016 for
Antarctica have been downloaded using MODIS package, then reprojected
and cropped to use them immediately in the examples.

Aspect and Digital Elevation Model of Antarctica obtained from National
Snow & Ice Center (<https://nsidc.org/data/nsidc-0082/versions/2>)\],
used as covariables after matching their extent and origin with
Aqua/Terra scenes.

``` r
library(modislst)
library(raster)
#> Loading required package: sp

## Add the target scene (aqua), and the reference scene (terra)
aqua <- raster(system.file("MODIS_data","aqua.tif", package="modislst"))
terra <- raster(system.file("MODIS_data","terra.tif", package="modislst"))


## Add covariables to improve the results 
dem <- raster(system.file("covariables","dem.tif", package="modislst"))
aspect <- raster(system.file("covariables","aspect.tif", package="modislst"))

## Plot the Aqua/Terra scene before the filling process
plot(aqua,main="Aqua Scene Before Filling - Target Scene")
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r
plot(terra, main="Terra Scene - Reference Scene")
```

<img src="man/figures/README-example-2.png" width="100%" />

### Fill the gaps in time

``` r
## Fill the gaps partially in time using moving window of size 47*47.
filledAqua <- fillingInTime(inputs = list(aqua,terra), m=47)

## Plot the filled Aqua scene 
plot(filledAqua, main="Filled Aqua Scene in time")
```

<img src="man/figures/README-reconstruct_in_time-1.png" width="100%" />

### Fill the gaps in space

``` r

## Fill the remaining gaps partially in space using.
filledAqua <- fillingInSpace(lstInput=filledAqua, covariable1=dem, covariable2=aspect)

## Plot the filled Aqua scene 
plot(filledAqua, main="Filled Aqua Scene")
```

<img src="man/figures/README-reconstruct_in_space-1.png" width="100%" />
