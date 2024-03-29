% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{modis_download}
\alias{modis_download}
\title{Download MODIS products}
\usage{
modis_download(
  dates = NULL,
  outdir = getwd(),
  tmpdir = tempdir(),
  username = NULL,
  password = NULL,
  messages = T,
  overwrite = F,
  overwrite_temp = F
)
}
\arguments{
\item{dates}{character vector of dates. Dates for which a modis image should
be downloaded. Needs to be of format "YYYY-mm-dd".}

\item{outdir}{character. Directory were the downloaded and processed files
should be stored to. This data will be used in \code{modis_classify()}}

\item{username}{character. Username of your EarthData account}

\item{password}{character. Password of your EarthData account}

\item{overwrite}{logical. If such a file already exists in the path, should
it be overwritten?}

\item{overwrite_temp}{logical. If the temporary HDF file already exists in
the tmpdir, should it be overwritten?}

\item{tempdir}{character. Directory were intermediary files should be stored
to (e.g. hdf files). By default these are stored to \code{tempdir()}.
However, you can change this to any other folder.}
}
\value{
Character vector of filename pointing to the downloaded files
}
\description{
This function is used to download and process modis files that can then be
classified into floodmaps using \code{modis_classify()}. Important: MODIS
bands 1-7 will be downloaded, yet only band 7 should be used in
\code{modis_classify()}. Note that the filenames simply correspond to the
aqcuisition date of the MODIS satellite imagery.
}
\examples{
\dontrun{
# Download files for two dates
files <- modis_download(
    dates     = c("2020-01-01", "2020-01-01")
  , outdir    = getwd()
  , tmpdir    = tempdir()
  , username  = "username"
  , password  = "password"
  , overwrite = F
)

# Check them
files
}
}
