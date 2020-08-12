################################################################################
#### DELETE
################################################################################
# Description: This script allows to download, process and classify MODIS
# imagery automatically.

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)
library(raster)
library(terra)
library(gdalUtils)
library(lubridate)
library(rgeos)
library(velox)
library(watermapr) # My own package -> see github

# Set the working directory
wd <- "/home/david/Schreibtisch"
setwd(wd)

# Make use of multiple cores
beginCluster()

# Specify an area of interest
aoi <- shapefile("/home/david/ownCloud/University/15. PhD/00_WildDogs/03_Data/01_RawData/MODIS/MCD43A4/02_Masks/OkavangoExtent")

# Specify your login credentials for EarthData
username <- "DoDx9"
password <- "EarthData99"

# Do you want to use dynamic watermasks? Note that this will drastically
# increase the time required to obtain floodmaps
dynamic <- TRUE

# Specify dates for which you want inundation maps
dates <- c(
    "2012-01-09"
  , "2012-01-17"
)

################################################################################
#### Level 1 Functions
################################################################################
#' Download MODIS products
#'
#' This function uses the data from the \code{modis_files()} function to
#' download the desired MODIS files.
#' @export
#' @param dataframe \code{data.frame}. Selection of the files that should be
#' downloaded (see \code{modis_files()})
#' @param path character. Path were the downloaded files should be stored to
#' @param username character. Username of your EarthData account
#' @param password character. Password of your EarthData account
#' @param overwrite logical. If such a file already exists in the path, should
#' it be overwritten?
#' @return Character vector of filename pointing to the downloaded files
modis_download <- function(
    date      = NULL
  , outdir    = getwd()
  , username  = NULL
  , password  = NULL
  , overwrite = F){

  # Identify all files that we need to download
  todownload <- lapply(1:length(date), function(x){

    # Show available files
    files <- modis_files(
        product     = "MCD43A4"
      , version     = "006"
      , start_date  = dates[x]
      , end_date    = dates[x]
      , aoi         = aoi
      , limit       = 100
    )

    # Subset to exact date (should be two files)
    files <- subset(files, as.character(date) == date[x])

    # Make sure there are exactly 2 files
    if (nrow(files) != 2){
      stop("There are more than two files!\n")
    }

    # Return dataframe of all files that need to be downloaded
    return(files)
  }) %>% do.call(rbind, .)

  # Loop through all files
  downloaded <- rep(NA, nrow(todownload))
  for (i in 1:nrow(todownload)){

    # And download the respective files
    downloaded[i] <- .modisDownload(
        dataframe = todownload[i, ]
      , path      = tempdir()
      , username  = username
      , password  = password
      , overwrite = overwrite
    )

    # Give a progress update
    cat(i, "out of", nrow(todownload), "tiles downloaded\n")
  }

  # Extract tiffs
  extracted <- lapply(1:length(downloaded), function(x){
    modis_extract(
        filepath  = downloaded[x]
      , outdir    = tempdir()
      , overwrite = overwrite
    )
  }) %>% do.call(c, .)

  # Group tiles that belong to the same date
  dat <- data.frame(Path = extracted, Date = modisDate(basename(extracted)))
  dat <- dat %>% group_by(date) %>% nest()

  # Stitch tiles
  stitched <- lapply(1:nrow(dat), function(x){
    modis_stitch(
        filepaths = dat$data[[x]]$Path
      , outdir    = outdir
      , outname   = paste0(dat$date[x], ".tif")
      , overwrite = overwrite
    )
  }) %>% do.call(c, .)

  # Reproject and crop stitched raster, then save
  final <- lapply(1:length(stitched), function(x){
    r <- rast(stitched[x])
    r <- project(r, y = CRS("+init=epsg:4326"), method = "near")
    r <- crop(r, aoi)
    names(r) <- paste0("Band_", 1:7)

    # Store it
    r <- writeRaster(
        r
      , filename  = stitched[x]
      , format    = "GTiff"
      , overwrite = TRUE
      , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
    )

    # Return the final directory
    return(stitched[x])
  }) %>% do.call(c, .)
  return(final)
}

#' Classify MODIS Image
#'
#' Function to classify a MODIS image into the binary categories water (valued
#' 1) and dryland (valued 0)
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param water \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent waters
#' @param dryland \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent dryland
#' @return \code{RasterLayer} of classified MODIS image
modis_classify <- function(x, water, dryland){

  # Check if the image is bimodal
  perc <- modis_percentiles(x, water = water, dryland = dryland)
  bimodal <- perc$WaterPercentiles[2] - 10 / 255 < perc$DrylandPercentiles[1]

  # If the image is not bimodal, we can return a message and skip the rest
  if (!bimodal){
    stop("Image not bimodal. Could not classify floodmap.")
  }

  # Make sure that the water and dryland masks have the same crs as the modis
  # layer
  water   <- spTransform(water, crs(x))
  dryland <- spTransform(dryland, crs(x))
  nowater <- spTransform(nowater, crs(x))

  # Extract the spectral values of band 7 below the dryland and water
  # polygons. To speed up the extraction we coerce the MODIS band to a velox
  # raster
  band7 <- velox(x)
  wat <- band7$extract(water) %>%
    do.call(rbind, .) %>%
    as.vector() %>%
    median(na.rm = TRUE)
  dry <- band7$extract(dryland) %>%
    do.call(rbind, .) %>%
    as.vector() %>%
    median(na.rm = TRUE)

  # Calculate the classification threshold
  mu <- wat + 0.3 * (dry - wat)

  # Predict/Classify the MODIS image into water (1) and dryland (0) using the
  # above calculated threshold (Clouds will become NA anyways)
  pred <- reclassify(x, c(0,mu,1, mu,1,0))

  # Get rid of areas that are most likely misclassified or only small ponds
  # that are irrelevant for our study. Set the raster values below this
  # polygon to 0
  pred <- raster::mask(pred, nowater, updatevalue = 0, inverse = TRUE)

  # Write the raster to a temporary file
  pred <- writeRaster(pred, tempfile())

  # Return the classified image
  return(pred)
}

################################################################################
#### Level 2 Functions
################################################################################
#' Show available MODIS products
#'
#' Function to show available products. This function allows you to determine
#' which MODIS products are available for download.
#' @export
#' @param product character indicating the desired products
#' @return \code{data.frame} containing all available products that fulfill the
#' search requirement. The column \code{short_name} can be used to retrieve
#' available files using the \code{modis_files()} function.
#' @examples
#' modis_products("MCD43A4")
modis_products <- function(product){

  # Retrieve available products from cmr
  prods <- read.csv("https://cmr.earthdata.nasa.gov/search/humanizers/report")

  # Remove undesired columns
  prods$original_value <- NULL
  prods$humanized_value <- NULL

  # Grep desired product
  if (!missing(product)){

    # If multiple products, paste them together into a single grep call
    if (length(product) > 1){
      product <- paste0(product, collapse = "|")
    }

    # Subset
    prods <- prods[grep(product, prods$short_name), ]
  }

  # Coerce everything to character
  prods <- dplyr::mutate_all(prods, as.character)

  # Return unique entries
  return(unique(prods))
}

#' Get Available Files of a MODIS Product
#'
#' Function to retrieve available files from a modis product. This function
#' allows you to find and select files that you want to download later using
#' \code{modis_download}.
#' @export
#' @param product character. Available products can be shown using
#' \code{modis_products}
#' @param version character. Indicates which version of the product should be
#' used. By default this is set to "006"
#' @param start_date character. Start date for the data requested formatted
#' yyyy-mm-dd
#' @param end_date character. End date for the data requested formatted
#' yyyy-mm-dd
#' @param aoi Any object from which an extent (using \code{bbox()} can be
#' derived)
#' @param limit positive integer. Number of results collected before the process
#' is aborted
#' @return \code{data.frame} containing some information about all files that
#' can be downloaded. The column "Online Access URLs" can be used in the
#' \code{modis_download()} function
#' @examples
#' files <- modis_files(product = "MCD43A4", version = "006", start_date =
#' "2020-02-01", end_date = "2020-02-04", aoi = '21.75,-20.65,24.3,-18.15',
#' limit = 100)
modis_files <-  function(
    product     = "MCD43A4"
  , version     = "006"
  , start_date  = NULL
  , end_date    = NULL
  , aoi         = NULL
  , limit       = 100){

  # Make sure the product exists
  dd <- modis_products(product)
  dd <- dd[dd$short_name == product & dd$version == version, ]

  # In case there is no such data, interrupt. In case there are multiple sources
  # of data, continue and use first
  if (nrow(dd) < 1){
    stop("Requested product not available")
  } else if (nrow(dd) > 1){
    warning("Multiple sources available, using first")
    print(dd)
    dd <- dd[1, ]
  }

  # Specify date suffix
  datesuffix <- "T00:00:00Z"

  # Convert dates to proper dates
  start_date  <- as.Date(start_date)
  end_date    <- as.Date(end_date)
  temporal    <- paste0(start_date, datesuffix, ",", end_date, datesuffix)

  # Put URLs together
  cmr_host <- "https://cmr.earthdata.nasa.gov"
  url <- file.path(cmr_host, "search/granules")

  # Retrieve extent
  ext <- .getExtent(aoi)

  # Put query parameters together
  params <- list(
      short_name    = product
    ,	temporal      = temporal
    , downloadable  = "true"
    , bounding_box  = ext
  )

  # Run query
  results <- .getSearchResults(url = url, limit = limit, args = params)

  # Identify modis dates
  results <- cbind(results, modisDate(results$`Producer Granule ID`))

  # Return the results
  return(results)
}

#' Extract Modis Bands and Convert HDF to Raster
#'
#' Function to convert the downloaded hdf files to proper .tif rasters
#' @export
#' @param filepath character. Filepath pointing to the downloaded hdf files.
#' Note that the filenames need to be unchanged and as downloaded from using
#' \code{modis_download()}
#' @param outdir character. Directory to which the converted file should be
#' stored. By default it is stored in the same folder as the hdf file.
#' @param removeHDF logical. Should the original hdf file be removed after
#' conversion?
#' @return filepath to the converted file
modis_extract <- function(
    filepath  = NULL
  , outdir    = dirname(filepath)
  , removeHDF = F
  , overwrite = F){

  # Extract the date from the modis file
  date <- as.vector(as.matrix(modisDate(filepath)))

  # Get all bands
  bands <- get_subdatasets(filepath)

  # Select the bands we want to keep
  bands <- bands[grep("Nadir.*Band[1-7]", bands)]

  # Create new filenames under which the bands will be stored
  names <- paste0(dirname(filepath), "/", basename(bands), ".tif")

  # Create output filename of the merged file
  filename <- paste0(outdir, "/", basename(filepath), ".tif")

  # Make sure the file does not already exist
  if (file.exists(filename) & !overwrite){
    warning(paste0("file ", filename, " already exists"))
    return(filename)
  }

  # Convert the selected bands to tifs using the new names
  for (j in 1:length(bands)){
    gdal_translate(bands[j], dst_dataset = names[j])
  }

  # Stack them
  r <- rast(names)

  # Assign date as band names
  names(r) <- paste0("Band_", 1:7)

  # Store the stacked file
  writeRaster(
      r
    , filename  = filename
    , format    = "GTiff"
    , overwrite = TRUE
    , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
  )

  # Remove the seperate bands
  file.remove(names)

  # In case the original HDF should be removed, do so
  if (removeHDF){
    file.remove(filepath)
  }

  # Return the location of the file
  return(filename)
}

#' Stitch Modis Tiles
#'
#' Function that allows you to stitch modis tiles together
#' @export
#' @param filepaths character. Filepaths pointing to the tiff files as prepared
#' using the \code{modis_extract} function. The selected files will be stitched
#' @param outname character. Filepath and filename to which the stitched tiff
#' should be stored
#' @return character pointing to the stitched tiff file
modis_stitch <- function(filepaths, outdir = getwd(), outname = NULL, overwrite = F){

  # Create output name
  filename <- paste0(outdir, "/", outname)
  if (file.exists(filename) & !overwrite){
    warning(paste0("file ", filename, " already exists"))
    return(filename)
  }

  # Create a virtual raster
  name <- tempfile(fileext = ".vrt")
  gdalbuildvrt(gdalfile = filepaths, output.vrt = name)

  # Coerce virtual raster to a true raster
  gdal_translate(
      src_dataset   = name
    , dst_dataset   = filename
    , output_Raster = TRUE
    , options       = c("BIGTIFFS=YES")
  )

  # Return the filepath to the stitched file
  return(filename)
}

#' Retrieve Date from MODIS Filenames
#'
#' Function that allows you to retrieve the date from a MODIS file
#' @export
#' @param filename character. Filename(s) of the files for which a date should
#' be extracted
#' @return character. Dates as derived from the MODIS filename
modisDate <- function(filename){
  ff <- basename(filename)
  dot <- sapply(strsplit(ff, "\\."), '[', 2)
  dates <- gsub("[aA-zZ]", "", dot)
  dates <- substr(basename(filename), 10, 16)
  dates <- .dateFromYearDoy(dates)
  data.frame(date = dates, stringsAsFactors = FALSE)
}

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
modis_watermask <- function(
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

#' Compare Spectral Signatures of MODIS Maps
#'
#' Function to compare spectral signatures of MODIS maps below water and dryland
#' polygons
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param water \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent waters
#' @param dryland \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent dryland
#' @return Plot of spatial reflectance values below the two polygons
modis_specs <- function(x = NULL, water = NULL, dryland = NULL){

    # Make sure that the water and dryland masks have the same crs as the modis
    # layer
    water   <- spTransform(water, crs(x))
    dryland <- spTransform(dryland, crs(x))

    # Extract the spectral values of band 7 below the dryland and water polygons
    # To speed up the extraction we coerce the modis band to a velox raster
    band7 <- velox(x)

    # Extract the values and coerce the output to a dataframe
    wat <- band7$extract(water) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    dry <- band7$extract(dryland) %>%
      do.call(rbind, .) %>%
      as.data.frame()

    # Prepare a column that indicates the land cover class
    wat$Class <- "Water"
    dry$Class <- "Dryland"

    # Bind the extracted values together
    specs <- rbind(wat, dry)

    # Plot the two densities for the spectral signatures of each value
    ggplot(specs, aes(V1, fill = Class)) + geom_density(alpha = 0.2)
}

#' Calculate Reflectance Percentiles
#'
#' Function that allows you to calculate percentiles in water and dryland
#' reflectance values
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param water \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent waters
#' @param dryland \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent dryland
#' @return dataframe of percentiles
modis_percentiles <- function(x, water, dryland){

  # Make sure that the water and dryland masks have the same crs as the modis
  # layer
  water   <- spTransform(water, crs(x))
  dryland <- spTransform(dryland, crs(x))

  # Extract the spectral values of band 7 below the dryland and water polygons.
  # To speed up the extraction we coerce the modis band to a velox raster
  band7 <- velox(x)

  # Extract the values and coerce the output to a dataframe
  wat <- band7$extract(water) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  dry <- band7$extract(dryland) %>%
    do.call(rbind, .) %>%
    as.data.frame()

  # Calculate the percentiles
  wat_perc <- quantile(wat[, 1], probs = c(0.01, 0.99), na.rm = TRUE)
  dry_perc <- quantile(dry[, 1], probs = c(0.01, 0.99), na.rm = TRUE)

  # Prepare a dataframe
  percs <- data.frame(WaterPercentiles = wat_perc, DrylandPercentiles = dry_perc)

  # Define output
  return(percs)
}

################################################################################
#### Helper Functions
################################################################################
# Helper function to download a single MODIS file
.modisDownload <- function(
    dataframe = NULL
  , path      = getwd()
  , username  = NULL
  , password  = NULL
  , overwrite = F){

  # Extract url
  url <- dataframe$`Online Access URLs`

  # Specify output directory
  filename <- file.path(path, basename(url))

  # We only download the file if it doesn't exist yet, or if we want to
  # overwrite it anyways
  if (!file.exists(filename) | overwrite){

    # Download file
    file <- httr::GET(
        url
      , httr::authenticate(username, password)
      , httr::progress()
      , httr::write_disk(filename, overwrite = overwrite)
    )
  } else {
    warning(paste0("file ", filename, " already exists, skipping download"))
  }
  return(filename)
}

# Helper function to prepare urls using query params
.getSearchResults <- function(
    url     = NULL
  , limit   = NULL
  , args    = NULL){

  # We start with the first page
  page_num <- 1

  # Initiate an empty results vector
  results <- NULL

  # We continue to read data until the predefined limit is reached
  while (length(results) < limit){

    # Read data using api
		response <- httr::GET(
        url   = url
      , query = c(args, page_num = page_num)
      , httr::add_headers(Accept = "text/csv")
    )

    # Check for a valid response
    httr::stop_for_status(response)

    # Make sure the response is of format "text/csv"
    if (httr::http_type(response) == "text/csv"){

      # Read content as dataframe
			data <- utils::read.csv(
          text             = httr::content(response, as = "text")
        , check.names      = FALSE
        , stringsAsFactors = FALSE
      )

      # Verify that the URL column is not empty
      catcher <- tryCatch(
          urls <- data[["Online Access URLs"]]
        , error = function(e){e}
      )

      # Make sure there is no error
      if (!inherits(catcher, "error")){
        if (length(urls) == 0){
          break
        }

        # Append the full table of results
        results <- rbind(results, data)
        page_num <- page_num + 1

      # In case there is an error
      } else {

        # We break the loop
        break
      }

    # If the response is not of format "text/csv"
    } else {

      # We break the loop
      break
    }
  }
  return(results)
}

# Helper function to retrieve the extent of a spatial object
.getExtent <- function(aoi){
  if (!is.vector(aoi)){
    b <- as.vector(bbox(aoi))
  } else {
    b <- as.vector(t(matrix(aoi, ncol = 2)))
  }
  paste(b, collapse = ",")
}

# Function to extract date from modis filename
.dateFromYearDoy <- function(x){
  year  <- as.integer(substr(x, 1, 4))
  doy   <- as.integer(substr(x, 5, 8))
  return(as.Date(doy, origin = paste(year - 1, "-12-31", sep = '')))
}
