
# floodmapr <img src="man/figures/floodmapr.png" align="right" width="150" height="150"/>

`floodmapr` is an R-package that allows you to download and classify
[MODIS MCD43A4](https://lpdaac.usgs.gov/products/mcd43a4v006/) satellite
imagery into binary maps of dryland and water cover. To be able to
download data, you need to have an
[EarthData](https://earthdata.nasa.gov/) account (free). The
classification algorithm is based on the publication of Wolski et al.,
2017 and currently only applicable for the extent of the Okavango Delta
(see [Okavango Research
Institute](http://www.okavangodata.ub.bw/ori/monitoring/flood_maps/#)).

## Installation

To install `floodmapr` you need to have the package `devtools`
installed. You can simply do:

``` r
library(devtools)
install_github("DavidDHofmann/floodmapr")
```

## Workflow

The `floodmapr` package follows a three step process using the
functions, `modis_download()`, `modis_load()`, and `modis_classify()`.

1.  `modis_download()` Download MCD43A4 satellite imagery for desired
    dates for the extent of the Okavango Delta. The function
    automatically downloads, extracts, stitches and reprojects (EPSG
    4326) corresponding modis tiles and stores them into a multiband
    raster. The output files are named according to the input date
    (e.g. *2020-01-01.tif*). Note: while only band 7 will be used in
    the later process, bands 1-7 are stored for completeness.

2.  `modis_load()` Once you have downloaded the satellite imagery, load
    it into R using `modis_load()`. This will only load band 7 into your
    current R session.

3.  `modis_classify()` Classify the satellite image into a binary map of
    water (valued 1) and dryland (valued 0). You can then use the
    resulting raster and e.g. `plot()` it.

In addition to the above described main functions there are some
additional helper functions that allow you to do some preliminary
analysis on MODIS MCD43A4 band 7.

  - `modis_bimodal()` allows you to check if modis band 7 is bimodal.
    This is required to be able to classify the image. Bimodality can
    fail during some months of the year, when surface reflectance of
    water and dryland become more alike.

  - `modis_specs()` plots the densities of water- and
    dryland-reflectance values. Again, this can be used to visualize how
    bimodal the distribution of these values is.

  - `modis_precentiles()` percentile values of water- and
    dryland-reflecatances.

Finally, there is a function called `modis_watermask()`. This function
can be used to create a dynamic watermask that can ultimately be fed
into `modis_classify()`. In contrast to the static watermask used by
`floodmapr` in the background, the dynamic watermask uses previously
classified foodmaps to determine areas that were frequently flooded in
the past. This way, the algorithm is able to learn about potential
long-term changes in the system.

## Example

Here is a complete example of the above outlined workflow

``` r
library(floodmapr)

# Download modis data
downloaded <- modis_download(
    dates     = c("2020-01-01", "2020-01-02")
  , outdir    = getwd()
  , tmpdir    = getwd()
  , username  = "username"
  , password  = "password"
  , overwrite = F
)

# Load data
loaded <- modis_load(downloaded[1])

# Do some checks
modis_specs(loaded)
modis_bimodal(loaded)
modis_percentiles(loaded)

# Classify data
classified <- modis_classify(loaded, ignore.bimodality = T)

# Visualize classified image
plot(classified)
```

## References

*Wolski, Piotr, Mike Murray-Hudson, Kgalalelo Thito, und Lin Cassidy.
„Keeping It Simple: Monitoring Flood Extent in Large Data-Poor
Wetlands Using MODIS SWIR Data“. International Journal of Applied Earth
Observation and Geoinformation 57 (Mai 2017): 224–34.
<https://doi.org/10.1016/j.jag.2017.01.005>.*
