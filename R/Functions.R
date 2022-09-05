################################################################################
#### Load Dependencies
################################################################################
#' @import raster terra ggplot2 sp
#' @importFrom dplyr group_by mutate_all
#' @importFrom magrittr %>%
#' @importFrom tidyr nest
#' @importFrom lubridate years ymd
#' @importFrom rgeos gBuffer
#' @importFrom RGISTools modSearch
#' @importFrom gdalUtils gdal_translate gdalbuildvrt
NULL

################################################################################
#### Level 1 Functions
################################################################################
#' Download MODIS products
#'
#' This function is used to download and process modis files that can then be
#' classified into floodmaps using \code{modis_classify()}. Important: MODIS
#' bands 1-7 will be downloaded, yet only band 7 should be used in
#' \code{modis_classify()}. Note that the filenames simply correspond to the
#' aqcuisition date of the MODIS satellite imagery.
#' @export
#' @param dates character vector of dates. Dates for which a modis image should
#' be downloaded. Needs to be of format "YYYY-mm-dd".
#' @param outdir character. Directory were the downloaded and processed files
#' should be stored to. This data will be used in \code{modis_classify()}
#' @param tempdir character. Directory were intermediary files should be stored
#' to (e.g. hdf files). By default these are stored to \code{tempdir()}.
#' However, you can change this to any other folder.
#' @param username character. Username of your EarthData account
#' @param password character. Password of your EarthData account
#' @param overwrite logical. If such a file already exists in the path, should
#' it be overwritten?
#' @param overwrite_temp logical. If the temporary HDF file already exists in
#' the tmpdir, should it be overwritten?
#' @return Character vector of filename pointing to the downloaded files
#' @examples
#' \dontrun{
#' # Download files for two dates
#' files <- modis_download(
#'     dates     = c("2020-01-01", "2020-01-01")
#'   , outdir    = getwd()
#'   , tmpdir    = tempdir()
#'   , username  = "username"
#'   , password  = "password"
#'   , overwrite = F
#' )
#'
#' # Check them
#' files
#'}
modis_download <- function(
    dates          = NULL
  , outdir         = getwd()
  , tmpdir         = tempdir()
  , username       = NULL
  , password       = NULL
  , messages       = T
  , overwrite      = F
  , overwrite_temp = F
  ) {

  # library(terra)
  # library(gdalUtils)
  # library(tidyverse)
  # library(rgdal)
  # library(rgeos)
  # library(lubridate)
  # library(RGISTools)
  # dates <- c("2020-01-01")
  # outdir <- "/home/david/Schreibtisch"
  # tmpdir <- "/home/david/Schreibtisch"
  # username <- "DoDx9"
  # password <- "EarthData99"
  # overwrite <- T
  # overwrite_temp <- T
  # messages <- T
  # load("/home/david/ownCloud/University/15. PhD/General/R-Packages/floodmapr/R/sysdata.rda")

  # Error messsages
  if (missing(dates)){stop("Provide dates")}
  if (missing(username)){stop("Provide username to access earth data")}
  if (missing(password)){stop("Provide password to access earth data")}
  if (!dir.exists(outdir)){stop("Specified output directory does not exist")}
  if (!dir.exists(tmpdir)){stop("Specified temporary directory does not exist")}

  # Parse dates
  dates <- ymd(dates)

  # Retrieve area of interest
  aoi <- masks_polygons[masks_polygons$Description == "aoi", ]

  # Make sure there are no duplicates
  if (length(dates) != length(unique(dates))){
    warning("Some dates are duplicated. Using only unique dates")
    dates <- unique(dates)
  }

  # Identify all files that we need to download. For each date there should only
  # be two (two tiles). Since the search also returns files slightly before or
  # after the desired date, we need to subset.
  todownload <- lapply(seq_along(dates), function(x) {
    files <- .modisFiles(
        product     = "MCD43A4"
      , version     = "006"
      , start_date  = dates[x]
      , end_date    = dates[x]
      , aoi         = aoi
    )
    if (nrow(files) != 2){
      stop("There are more than two files!\n")
    }
    return(files)
  }) %>% do.call(rbind, .)

  # Check if some of the files already exist
  files <- file.path(tmpdir, basename(todownload$URL))

  # Download all selected files to a temporary directory
  downloaded <- rep(NA, nrow(todownload))
  for (i in 1:nrow(todownload)) {
    if (messages) {
      cat("Downloading tiles:", i, "out of", nrow(todownload), "\n")
    }
    if (file.exists(files[i]) & !overwrite_temp) {
      if (messages) {
        cat("File", files[i], "exists and is not overwritten\n")
      }
      downloaded[i] <- files[i]
    } else {
      downloaded[i] <- .modisDownload(
          dataframe = todownload[i, ]
        , path      = tmpdir
        , username  = username
        , password  = password
        , overwrite = T
      )
    }
  }

  # Extract tiffs
  if (messages) {
    cat("All tiles downloaded. Extracting tiffs from hdf files now...\n")
  }
  extracted <- lapply(1:length(downloaded), function(x) {
    .modisExtract(
        filepath  = downloaded[x]
      , outdir    = tmpdir
      , overwrite = overwrite_temp
      , removeHDF = F
    )
  }) %>% do.call(c, .)

  # There are two tiles per date. We need to stitch them. Thus, create a tibble
  # that shows which files need to be combined
  if (messages) {
    cat("All hdf files extracted. Stitching tiles now...\n")
  }
  dat <- data.frame(
      Path = extracted
    , Date = modis_date(basename(extracted)) %>% as.matrix() %>% as.vector()
  ) %>% mutate_all(., as.character)
  dat <- dat %>% group_by(Date) %>% nest()
  stitched <- lapply(1:nrow(dat), function(x) {
    suppressWarnings(
      .modisStitch(
          filepaths = dat$data[[x]]$Path
        , outdir    = outdir
        , outname   = paste0(dat$Date[x], ".tif")
        , overwrite = overwrite
      )
    )
  }) %>% do.call(c, .)

  # Remove unstitched files
  file.remove(extracted)

  # Reproject and crop stitched raster, then save. Note that we will use ngb for
  # resampling. The reason is that bilinear interpolation would potentially
  # contaminate or be contaminated by NA values.
  if (messages) {
    cat("All tiles stitched. Reprojecting and cropping final tiles now...\n")
  }
  final <- lapply(1:length(stitched), function(x){
    r <- rast(stitched[x])
    r <- suppressWarnings(crop(r, spTransform(aoi, crs(r)), snap = "out"))
    r <- terra::project(r, "+proj=longlat +datum=WGS84 +no_defs", method = "near")
    r <- crop(r, aoi, snap = "out")
    names(r) <- paste0("Band_", 1:7)
    r <- terra::writeRaster(
        r
      , filename  = stitched[x]
      , filetype  = "GTiff"
      , overwrite = TRUE
      , gdal      = c("INTERLEAVE = BAND", "COMPRESS = LZW")
    )
    return(stitched[x])
  }) %>% do.call(c, .)

  # Return the location of the final file
  if (messages) {
    cat("Finished!\n")
  }
  return(final)
}

#' Load Downloaded MODIS Data
#'
#' Function to load downloaded MODIS MCD43A4 band 7 (see
#' \code{modis_download()}). This function simply wraps the command
#' \code{raster(filepath, band = 7)}.
#' @export
#' @param filepath character. Filepath pointing to the downloaded modis file
#' @return \code{RasterLayer} containing only MODIS MCD43A4 band 7
#' @examples
#' \dontrun{
#' # Download files for two dates
#' files <- modis_download(
#'     dates     = c("2020-01-01", "2020-01-01")
#'   , outdir    = getwd()
#'   , tmpdir    = tempdir()
#'   , username  = "username"
#'   , password  = "password"
#'   , overwrite = F
#' )
#'
#' # Load one of them and show it
#' modis <- modis_load(files[1])
#' show(modis)
#'}
modis_load <- function(filepath = NULL) {

  # Make sure only one file is given
  if (length(filepath) > 1){
    stop("Can't load multiple files at once. Provide a single file only")
  }

  # Load it
  return(rast(filepath)[[7]])
}

#' Classify MODIS Image
#'
#' Function to classify MODIS MCD43A4 band 7 into the binary categories water
#' and dryland. Classification is based on the algorithm described in Wolski et
#' al., 2017 and requires that reflectance values of water- and dryland are
#' sufficiently distinct. The final map binarily depicts water = 1 and dryland =
#' 0.
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param watermask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the water-polygon that is used to extract reflectances of water
#' on MODIS band 7. By default the masks from Wolski et al., 2017 are used. See
#' also \code{modis_watermask()} in case you want to create dynamic watermasks.
#' @param drymask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the dryland-polygon that is used to extract reflectances of
#' water on MODIS band 7. By default the masks from Wolski et al., 2017 are
#' used.
#' @param ignore.bimodality logical. Should issues with bimodality be ignored,
#' i.e. the bimodality check be skipped? This can lead to biased
#' classifications but may help in detecting issues.
#' @return \code{RasterLayer} of classified MODIS image. water is valued 1,
#' dryland valued 0. If there are clouds, they are masked as NA.
#' @examples
#' \dontrun{
#' # Download files for two dates
#' files <- modis_download(
#'     dates     = c("2020-01-01", "2020-01-01")
#'   , outdir    = getwd()
#'   , tmpdir    = tempdir()
#'   , username  = "username"
#'   , password  = "password"
#'   , overwrite = F
#' )
#'
#' # Load one of them
#' modis <- modis_load(files[1])
#'
#' # Classify it
#' classified <- modis_classify(modis, ignore.bimodality = T)
#'
#' # Visualize
#' plot(classified)
#'}
modis_classify <- function(
    x                 = NULL
  , watermask         = NULL
  , drymask           = NULL
  , ignore.bimodality = F) {

  # Retrieve water, nowater and dryland masks
  water   <- vect(masks_polygons[masks_polygons$Description == "water", ])
  dryland <- vect(masks_polygons[masks_polygons$Description == "dryland", ])
  nowater <- vect(masks_polygons[masks_polygons$Description == "nowater", ])

  # In case a watermask or drymask is provided, use them
  if (!is.null(watermask)){water <- watermask}
  if (!is.null(drymask)){dryland <- drymask}

  # Check if the image is bimodal
  if (!ignore.bimodality) {
    bimodal <- modis_bimodal(
        x         = x
      , watermask = water
      , drymask   = dryland
    )

    # If the image is not bimodal, we can return a message and skip the rest
    if (!bimodal){
      stop("Image not bimodal. Could not classify floodmap.")
    }
  }

  # Extract the spectral values of band 7 below the dryland and water
  # polygons
  wat <- mask(x, water) %>% as.vector() %>% na.omit() %>% median()
  dry <- mask(x, dryland) %>% as.vector() %>% na.omit() %>% median()

  # Calculate the classification threshold
  mu <- wat + 0.3 * (dry - wat)

  # Predict/Classify the MODIS image into dryland (0), water (1), and clouds
  # (2), using the above calculated threshold
  # pred <- reclassify(x, c(0,mu,1, mu,1,0))
  rcl <- data.frame(from = c(NA, 0, mu), to = c(NA, mu, 1), new = c(2, 1, 0))
  pred <- classify(x, rcl)

  # Get rid of areas that are always dry
  pred <- mask(pred, nowater, updatevalue = 0, inverse = TRUE)

  # Write the raster to a temporary file
  pred <- terra::writeRaster(pred, tempfile(fileext = ".tif"))

  # Return the classified image
  return(pred)
}

################################################################################
#### Level 2 Functions
################################################################################
#' Check for Bimodality in MODIS Data
#'
#' Function to check for bimodality in MODIS MCD43A4 band 7. Bimodality is said
#' to be achieved if peaks of water- and dryland reflectances are sufficiently
#' distinct. For details check Wolksi et al. 2017.
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param watermask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the water-polygon that is used to extract reflectances of water
#' on MODIS band 7. By default the masks from Wolski et al., 2017 are used. See
#' also \code{modis_watermask()}.
#' @param drymask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the dryland-polygon that is used to extract reflectances of
#' water on MODIS band 7. By default the masks from Wolski et al., 2017 are
#' used.
#' @return logical indicating if the image is bimodal (TRUE) or not bimodal
#' (FALSE)
#' @examples
#' \dontrun{
#' # Download file
#' file <- modis_download(
#'     dates     = "2020-01-01"
#'   , outdir    = getwd()
#'   , tmpdir    = tempdir()
#'   , username  = "username"
#'   , password  = "password"
#'   , overwrite = F
#' )
#'
#' # Load it
#' modis <- modis_load(file)
#'
#' # Check for bimodality
#' modis_bimodal(modis)
#'}
modis_bimodal <- function(
    x                 = NULL
  , watermask         = NULL
  , drymask           = NULL) {

  # Retrieve water and dryland masks
  water   <- vect(masks_polygons[masks_polygons$Description == "water", ])
  dryland <- vect(masks_polygons[masks_polygons$Description == "dryland", ])

  # In case a watermask or drymask is provided, use them
  if (!is.null(watermask)) {water <- watermask}
  if (!is.null(drymask)) {dryland <- drymask}

  # Check if the image is bimodal
  perc <- modis_percentiles(x, watermask = water, drymask = dryland)
  bimodal <- perc$WaterPercentiles[2] - 10 / 255 < perc$DrylandPercentiles[1]

  # Return answer
  return(bimodal)
}

#' Retrieve Date from MODIS Filenames
#'
#' Function to retrieve the date from a MODIS hdf file. Note that this function
#' only works with the original MODIS hdf files.
#' @export
#' @param filename character. Filename(s) of the files for which a date should
#' be extracted
#' @return character. Dates as derived from the MODIS filename
modis_date <- function(filename = NULL) {
  ff <- basename(filename)
  dot <- sapply(strsplit(ff, "\\."), '[', 2)
  dates <- gsub("[aA-zZ]", "", dot)
  dates <- substr(basename(filename), 10, 16)
  dates <- .dateFromYearDoy(dates)
  data.frame(date = dates, stringsAsFactors = FALSE)
}

#' Create Dynamic Watermask
#'
#' Function to create a dynamic watermask for a point in time based on previous
#' watermaps. The function uses floodmaps from previous dates to determine areas
#' that were covered by water in at least x% of the time. This can help to deal
#' with issues of non-bimodality.
#' @export
#' @param date date. Date for which a watermask should be calculated
#' @param floodmaps character vector. Filenames of all floodmaps that already
#' exist
#' @param filedates dates vector. Dates which the above floodmaps represent
#' @param years numeric. Number of years that should be considered to create a
#' watermask
#' @param threshold numeric. Needs to be between 0 and 1. How often (relative
#' frequency) a pixel needs to be inundated in the past x years to be considered
#' in the watermask
#' @return \code{SpatialPolygons} of the waterask
modis_watermask <- function(
    date      = NULL
  , filenames = NULL
  , filedates = NULL
  , years     = 5
  , threshold = 0.99){

  # filenames <- dir("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/00_Floodmaps", full.names = T)
  # filedates <- ymd(basename(filenames))
  # date <- ymd("2018-01-01")
  # years <- 5
  # threshold <- 0.99

  # Some error messages
  if (missing(filenames)){ stop("Provide filenames") }
  if (missing(filedates)){ stop("Provide filedates") }
  if (threshold < 0 | threshold > 1){ stop("Threshold needs to be between 0 and 1") }
  if (years < 0){ stop("Can't have negative years") }

  # Make naming nicer
  end_date <- date

  # Subtract 5 years, to get the first date we would include to calculate the mask
  start_date <- end_date - years(5)

  # Identify all possible dates between start and end dates for which we would
  # include maps to calculate the mask
  period <- seq(start_date, end_date, "days")

  # Keep only those filenames which are within the period of interest
  filenames <- filenames[filedates %in% period]
  if (length(filenames) == 0){
    stop("No files available for the desired years")
  }

  # Load the files into a stack
  formask <- rast(filenames)

  # Reclassify the stack so that water becomes 1, dryland and clouds 0
  rcl <- data.frame(old = c(0, 1, 2), new = c(0, 1, 0))
  formask <- classify(formask, rcl)

  # Sum the layers
  sum <- sum(formask)

  # Identify areas where there was water 99% of the time
  areas <- sum > threshold * nlyr(formask)

  # Coerce SpatRaster to RasterLayer
  areas <- raster(areas)

  # Polygonize
  wetmask <- rasterToPolygons(
      areas
    , fun       = function(x){ x == 1 }
    , dissolve  = TRUE
  )

  # Apply a small negative buffer to avoid errors due to pixel size
  wetmask <- suppressWarnings(gBuffer(wetmask, width = -1 / 111 * 0.25))

  # Return the final watermask
  return(vect(wetmask))
}

#' Compare Spectral Signatures of MODIS Maps
#'
#' Function to compare spectral signatures of MODIS maps below water and dryland
#' polygons
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param watermask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the water-polygon that is used to extract reflectances of water
#' on MODIS band 7. By default the masks from Wolski et al., 2017 are used. See
#' also \code{modis_watermask()} in case you want to create dynamic watermasks.
#' @param drymask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the dryland-polygon that is used to extract reflectances of
#' water on MODIS band 7. By default the masks from Wolski et al., 2017 are
#' used.
#' @return Plot of spatial reflectance values below the two polygons
#' @examples
#' \dontrun{
#' # Download file
#' file <- modis_download(
#'     dates     = "2020-01-01"
#'   , outdir    = getwd()
#'   , tmpdir    = tempdir()
#'   , username  = "username"
#'   , password  = "password"
#'   , overwrite = F
#' )
#'
#' # Load it
#' modis <- modis_load(file)
#'
#' # Check spectral reflectances
#' modis_specs(modis)
#'}
modis_specs <- function(
    x         = NULL
  , watermask = NULL
  , drymask   = NULL){

  # Retrieve water and dryland masks
  water   <- vect(masks_polygons[masks_polygons$Description == "water", ])
  dryland <- vect(masks_polygons[masks_polygons$Description == "dryland", ])

  # In case a watermask or drymask is provided, use them
  if (!is.null(watermask)){water <- watermask}
  if (!is.null(drymask)){dryland <- drymask}

  # Extract the spectral values of band 7 below the dryland and water polygons
  wat <- mask(x, water) %>% as.data.frame() %>% na.omit()
  dry <- mask(x, dryland) %>% as.data.frame() %>% na.omit()

  # Prepare a column that indicates the land cover class
  wat$Class <- "Water"
  dry$Class <- "Dryland"

  # Bind the extracted values together
  specs <- rbind(wat, dry)
  specs <- na.omit(specs)
  names(specs)[1] <- "Band_7"

  # Plot the two densities for the spectral signatures of each value
  ggplot(specs, aes(Band_7, fill = Class)) + geom_density(alpha = 0.2) +
    theme_minimal() +
    scale_fill_manual(values = c("darkgreen", "cornflowerblue")) +
    xlab("Reflectance") +
    ylab("Denstiy") +
    ggtitle("MODIS Terra Surface Reflectance")
}

#' Calculate Reflectance Percentiles
#'
#' Function to calculate percentiles of water and dryland reflectance values
#' @export
#' @param x \code{RasterLayer} of MODIS band 7
#' @param watermask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the water-polygon that is used to extract reflectances of water
#' on MODIS band 7. By default the masks from Wolski et al., 2017 are used. See
#' also \code{modis_watermask()} in case you want to create dynamic watermasks.
#' @param drymask \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame}
#' representing the dryland-polygon that is used to extract reflectances of
#' water on MODIS band 7. By default the masks from Wolski et al., 2017 are
#' used.
#' @return dataframe of percentiles
#' @examples
#' \dontrun{
#' # Download file
#' file <- modis_download(
#'     dates     = "2020-01-01"
#'   , outdir    = getwd()
#'   , tmpdir    = tempdir()
#'   , username  = "username"
#'   , password  = "password"
#'   , overwrite = F
#' )
#'
#' # Load it
#' modis <- modis_load(file)
#'
#' # Check percentiles
#' modis_percentiles(modis)
#'}
modis_percentiles <- function(
    x         = NULL
  , watermask = NULL
  , drymask   = NULL) {

  # Retrieve water and dryland masks
  water   <- vect(masks_polygons[masks_polygons$Description == "water", ])
  dryland <- vect(masks_polygons[masks_polygons$Description == "dryland", ])

  # In case a watermask or drymask is provided, use them
  if (!is.null(watermask)){water <- watermask}
  if (!is.null(drymask)){dryland <- drymask}

  # Extract the spectral values of band 7 below the dryland and water polygons
  wat <- mask(x, water) %>% as.vector() %>% na.omit()
  dry <- mask(x, dryland) %>% as.vector() %>% na.omit()

  # Calculate the percentiles
  wat_perc <- quantile(wat, probs = c(0.01, 0.99), na.rm = TRUE)
  dry_perc <- quantile(dry, probs = c(0.01, 0.99), na.rm = TRUE)

  # Prepare a dataframe
  percs <- data.frame(WaterPercentiles = wat_perc, DrylandPercentiles = dry_perc)

  # Define output
  return(percs)
}

################################################################################
#### Helper Functions
################################################################################
# Function to extract information from modis filenames
.modisInfo <- function(urlname){
  name <- basename(urlname)
  info <- strsplit(name, split = "\\.")
  info <- lapply(info, rbind)
  info <- do.call(rbind, info)
  info <- as.data.frame(info, stringsAsFactors = F)
  names(info) <- c("Product", "AcquisitionDate", "Tile", "Version", "ProductionDate", "Format")
  info$AcquisitionDate <- substr(info$AcquisitionDate, start = 2, stop = nchar(info$AcquisitionDate))
  info$AcquisitionDate <- as.Date(info$AcquisitionDate, "%Y%j")
  info$ProductionDate <- as.Date(as.POSIXct(info$ProductionDate, format = "%Y%j%H%M%S"))
  return(info)
}

# Function to find available modis files
.modisFiles <- function(product = "MCD43A4", version = "006", start_date, end_date, aoi) {

  # Get available products
  search_res <- modSearch(
      product    = product
    , collection = version
    , startDate  = start_date
    , endDate    = end_date
    , extent     = extent(aoi)
  )
  search_res <- unname(search_res$hdf)

  # Extract file info
  results <- cbind(URL = search_res, .modisInfo(search_res), stringsAsFactors = F)

  # Subset to dates of interest
  results <- subset(results, AcquisitionDate >= start_date & AcquisitionDate <= end_date)

  # Return dataframe of results
  return(results)
}

# Helper function to download a single MODIS file
.modisDownload <- function(dataframe, path = getwd(), username, password, overwrite = F) {

  # Extract url
  url <- dataframe$URL

  # Specify output directory
  filename <- file.path(path, basename(url))

  # We only download the file if it doesn't exist yet, or if we want to
  # overwrite it anyways
  if ((!file.exists(filename)) | overwrite){

      # # Download file
      # download.file(
      #     url      = url
      #   , destfile = filename
      #   , method   = "wget"
      #   , extra    = paste("--user", username, "--password", password)
      # )

      # Download file
      file <- httr::GET(
          url
        , httr::authenticate(username, password, type = "any")
        , httr::progress()
        , httr::write_disk(filename, overwrite = overwrite)
      )
    } else {
      warning(paste0("file ", filename, " already exists, skipping download"))
  }
  return(filename)
}

# Function to extract date from modis filename
.dateFromYearDoy <- function(x){
  year  <- as.integer(substr(x, 1, 4))
  doy   <- as.integer(substr(x, 5, 8))
  return(as.Date(doy, origin = paste(year - 1, "-12-31", sep = '')))
}

# Function to convert the downloaded hdf files to proper .tif rasters
.modisExtract <- function(
    filepath  = NULL
  , outdir    = dirname(filepath)
  , removeHDF = F
  , overwrite = F) {

  # Load the layers
  bands <- sds(filepath)

  # Keep only the bands of interest
  bands <- bands[[grep("Nadir.*Band[1-7]", names(bands))]]

  # Convert to terra raster
  bands <- rast(bands)

  # Assign date as band names
  names(bands) <- paste0("Band_", 1:7)

  # Create output filename of the new file
  filename <- paste0(file.path(outdir, basename(filepath)), ".tif")

  # Make sure the file does not already exist, then save
  if (file.exists(filename) & !overwrite) {
    cat(paste0("file", filename, "already exists and will not be overwritten...\n"))
  } else {
    terra::writeRaster(bands, filename  = filename, overwrite = TRUE)
  }

  # In case the original HDF should be removed, do so
  if (removeHDF) {
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
  , overwrite = F) {

  # Create output name
  filename <- file.path(outdir, outname)
  if (file.exists(filename) & !overwrite){
    cat("file", filename, "already exists and is not overwritten...\n")
  } else {

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
  }

  # Remove aux files
  remove <- paste0(filename, ".aux.xml")
  file.remove(remove)

  # Return the filepath to the stitched file
  return(filename)
}
