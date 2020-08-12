################################################################################
#### Load Dependencies
################################################################################
#' @import raster terra velox ggplot2 sp
#' @importFrom dplyr group_by mutate_all
#' @importFrom magrittr %>%
#' @importFrom tidyr nest
#' @importFrom lubridate years
#' @importFrom rgeos gBuffer
#' @importFrom gdalUtils get_subdatasets gdal_translate gdalbuildvrt
NULL

################################################################################
#### Level 1 Functions
################################################################################
#' Download MODIS products
#'
#' This function is used to download and process modis files that can then be
#' classified into floodmaps using \code{modis_classify}. Important: MODIS bands
#' 1-7 will be downloaded, yet only band 7 should be used in
#' \code{modis_classify}. Note that the filenames simply correspond to the
#' aqcuisition date of the MODIS satellite imagery.
#' @export
#' @param date character vector of dates. Dates for which a modis image should be
#' downloaded. Needs to be of format "YYYY-mm-dd".
#' @param outdir character. Directory were the processed files should be stored
#' to. This is the raw data that will be used for floodmapping
#' @param tempdir character. Directory were intermediary files should be stored.
#' By default this is determined by \code{tempdir()}. However, you can change
#' this to another folder.
#' @param username character. Username of your EarthData account
#' @param password character. Password of your EarthData account
#' @param overwrite logical. If such a file already exists in the path, should
#' it be overwritten?
#' @return Character vector of filename pointing to the downloaded files
modis_download <- function(
    dates     = NULL
  , outdir    = getwd()
  , tmpdir    = tempdir()
  , username  = NULL
  , password  = NULL
  , overwrite = F){

  # Some error messsages
  if (missing(dates)){stop("Provide dates")}
  if (missing(username)){stop("Provide username to access earth data")}
  if (missing(password)){stop("Provide password to access earth data")}

  # Retrieve area of interest
  aoi <- masks_polygons[masks_polygons$Description == "aoi", ]

  # Identify all files that we need to download. For each date there should only
  # be two (two tiles). Since the search also returns files slightly before or
  # after the desired date, we need to subset.
  todownload <- lapply(1:length(dates), function(x){
    files <- .modisFiles(
        product     = "MCD43A4"
      , version     = "006"
      , start_date  = dates[x]
      , end_date    = dates[x]
      , aoi         = aoi
      , limit       = 100
    ) %>% subset(as.character(date) == dates[x])
    if (nrow(files) != 2){
      stop("There are more than two files!\n")
    }
    return(files)
  }) %>% do.call(rbind, .)

  # Download all selected files to a temporary directory
  downloaded <- rep(NA, nrow(todownload))
  for (i in 1:nrow(todownload)){
    cat("Downloading tiles:", i, "out of", nrow(todownload), "\n")
    downloaded[i] <- .modisDownload(
        dataframe = todownload[i, ]
      , path      = tmpdir
      , username  = username
      , password  = password
      , overwrite = overwrite
    )
  }

  # Extract tiffs
  cat("All tiles downloaded. Extracting tiffs from hdf files now...\n")
  extracted <- lapply(1:length(downloaded), function(x){
    .modisExtract(
        filepath  = downloaded[x]
      , outdir    = tmpdir
      , overwrite = overwrite
      , removeHDF = F
    )
  }) %>% do.call(c, .)

  # There are two tiles per date. We need to stitch them. Thus, create a tibble
  # that shows which files need to be combined
  cat("All hdf files extracted. Stitching tiles now...\n")
  dat <- data.frame(Path = extracted, Date = modis_date(basename(extracted)))
  dat <- dat %>% group_by(date) %>% nest()
  stitched <- lapply(1:nrow(dat), function(x){
    .modisStitch(
        filepaths = dat$data[[x]]$Path
      , outdir    = outdir
      , outname   = paste0(dat$date[x], ".tif")
      , overwrite = overwrite
    )
  }) %>% do.call(c, .)

  # Reproject and crop stitched raster, then save
  cat("All tiles stitched. Reprojecting and cropping final tiles now...\n")
  final <- lapply(1:length(stitched), function(x){
    r <- rast(stitched[x])
    r <- project(r, y = CRS("+init=epsg:4326"), method = "near")
    r <- crop(r, aoi)
    names(r) <- paste0("Band_", 1:7)

    # Store it and return the directory
    r <- writeRaster(
        r
      , filename  = stitched[x]
      , format    = "GTiff"
      , overwrite = TRUE
      , options   = c("INTERLEAVE = BAND", "COMPRESS = LZW")
    )
    return(stitched[x])
  }) %>% do.call(c, .)
  cat("Finished!")
  return(final)
}

#' Classify MODIS Image
#'
#' Function to classify a MODIS image into the binary categories water and
#' dryland. Water = 1, Dryland = 0.
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param watermask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the water-polygon that is used to extract reflectances of water
#' on MODIS band 7.
#' @param drymask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the dryland-polygon that is used to extract reflectances of
#' water on MODIS band 7.
#' @param ignore.bimodality logical. Should the bimodality check be ignored?
#' @return \code{RasterLayer} of classified MODIS image
modis_classify <- function(
    x                 = NULL
  , watermask         = NULL
  , drymask           = NULL
  , ignore.bimodality = F){

  # Retrieve water, nowater and dryland masks
  water   <- masks_polygons[masks_polygons$Description == "water", ]
  dryland <- masks_polygons[masks_polygons$Description == "dryland", ]
  nowater <- masks_polygons[masks_polygons$Description == "nowater", ]

  # In case a watermask or drymask is provided, use them
  if (!is.null(watermask)){water <- watermask}
  if (!is.null(drymask)){dryland <- drymask}

  # Check if the image is bimodal
  if (!ignore.bimodality){
    perc <- modis_percentiles(x, water = water, dryland = dryland)
    bimodal <- perc$WaterPercentiles[2] - 10 / 255 < perc$DrylandPercentiles[1]

    # If the image is not bimodal, we can return a message and skip the rest
    if (!bimodal){
      stop("Image not bimodal. Could not classify floodmap.")
    }
  }

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
#' Retrieve Date from MODIS Filenames
#'
#' Function to retrieve the date from a MODIS file. Note that this function only
#' works with the original MODIS filenames
#' @export
#' @param filename character. Filename(s) of the files for which a date should
#' be extracted
#' @return character. Dates as derived from the MODIS filename
modis_date <- function(filename){
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

    # Retrieve water and dryland maps
    water <- masks_polygons[masks_polygons$Description == "water", ]
    dryland <- masks_polygons[masks_polygons$Description == "dryland", ]

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
    specs <- na.omit(specs)

    # Plot the two densities for the spectral signatures of each value
    ggplot(specs, aes(V1, fill = Class)) + geom_density(alpha = 0.2)
}

#' Calculate Reflectance Percentiles
#'
#' Function to calculate percentiles in water and dryland reflectance values
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param water \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent waters
#' @param dryland \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing permanent dryland
#' @return dataframe of percentiles
modis_percentiles <- function(x, water, dryland){

  # Retrieve water and dryland maps
  water <- masks_polygons[masks_polygons$Description == "water", ]
  dryland <- masks_polygons[masks_polygons$Description == "dryland", ]

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

# Function to show available products. This function allows you to determine
# which MODIS products are available for download.
.modisProducts <- function(product){

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
  prods <- mutate_all(prods, as.character)

  # Return unique entries
  return(unique(prods))
}

# Function to retrieve available files from a modis product. This function
# allows you to find and select files that you want to download.
.modisFiles <-  function(
    product     = "MCD43A4"
  , version     = "006"
  , start_date  = NULL
  , end_date    = NULL
  , aoi         = NULL
  , limit       = 100){

  # Make sure the product exists
  dd <- .modisProducts(product)
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
  results <- cbind(results, modis_date(results$`Producer Granule ID`))

  # Return the results
  return(results)
}

# Function to convert the downloaded hdf files to proper .tif rasters
.modisExtract <- function(
    filepath  = NULL
  , outdir    = dirname(filepath)
  , removeHDF = F
  , overwrite = F){

  # Extract the date from the modis file
  date <- as.vector(as.matrix(modis_date(filepath)))

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

# Function to stitch modis tiles together
.modisStitch <- function(
    filepaths = NULL
  , outdir    = getwd()
  , outname   = NULL
  , overwrite = F){

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
