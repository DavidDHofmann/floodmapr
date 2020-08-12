################################################################################
#### Function to Create Dynamic Watermask
################################################################################
#' Create Dynamic Watermasks
#'
#' Function to createa dynamic watermask based on previous floodmaps.
#' @export
#' @param date date. Date for which a watermask should be calculated
#' @param floodmaps character. Filenames of all floodmaps that already exist
#' @param filedates vector of dates. Dates which the above floodmaps represent
#' @param years numeric. Number of years that should be considered to create a
#' watermask
#' @param Threshold numeric between 0 and 1. How often (relative frequency) a
#' pixel needs to be inundated in the past x years to be considered in the
#' watermask
#' @return \code{SpatialPolygons} of the waterask
waterMask <- function(
    date      = NULL
  , filenames = NULL
  , filedates = NULL
  , years     = 5
  , threshold = 0.99){

  # Make naming nicer
  end_date <- date

  # Subtract 5 years, to get the first date we would include to calculate the mask
  start_date <- end_date - years(5)

  # Identify all possible dates between start and end dates for which we would
  # include maps to calculate the mask
  period <- seq(start_date, end_date, "days")

  # Keep only those filenames which are within the period of interest
  filenames <- filenames[filedates %in% period]

  # Load the files into a stack
  formask <- rast(filenames, bands = 1)

  # Reclassify the stack so that water becomes 1, dryland and clouds 0
  rcl <- data.frame(old = c(0, 127, 255), new = c(1, 0, 0))
  formask <- classify(formask, rcl)

  # Sum the layers
  sum <- sum(formask)

  # Identify areas where there was water 99% of the time
  areas <- sum > threshold * nlyr(formask)
  areas <- raster(areas)

  # Polygonize
  wetmask <- rasterToPolygons(areas
    , fun = function(x){x == 1}
    , dissolve = TRUE
  )

  # Apply a small negative buffer to avoid errors due to pixel size
  wetmask <- suppressWarnings(gBuffer(wetmask, width = -1/111*0.25))

  # Return the final watermask
  return(wetmask)
}
