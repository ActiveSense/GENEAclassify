% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataImport.R
\name{dataImport}
\alias{dataImport}
\title{Data import function}
\usage{
dataImport(
  binfile,
  downsample = 100,
  start = NULL,
  end = NULL,
  Use.Timestamps = FALSE,
  radians = FALSE,
  keep_raw_data = TRUE,
  mmap.load = (.Machine$sizeof.pointer >= 8),
  ...
)
}
\arguments{
\item{binfile}{File path to binary data to be segmented.}

\item{downsample}{Rate to downsample the data, defaults to every 100th observation. For no downsampling set NULL.}

\item{start}{Where to start reading observations.}

\item{end}{Where to end reading observations.}

\item{Use.Timestamps}{To use timestamps as the startand end time values this has to be set to TRUE. (Default FALSE)}

\item{radians}{calculate degrees rotation in radians.}

\item{keep_raw_data}{Keep the raw data for calculating steps using stepcounter.}

\item{mmap.load}{Default is (.Machine$sizeof.pointer >= 8). see mmap for more details}

\item{...}{additional arguments passed through.}
}
\value{
Returns a list containing a matrix of the data including the x, y and z axis data, vectors of the up down (elevation) 
and degrees (rotation), a vector of time stamps, a vector of vector magnitudes and the serial number of the device.
}
\description{
Loads the data into R and creates format required for segmentation.
}
\details{
Reads in the binary data file and extracts the information required for the segmentation procedure.
}
\examples{
##   segData = dataImport(bindata = file.path(system.file(package = "GENEAread"),
##                                                         "binfile",
##                                                         "TESTfile.bin"))
## 
## segData1 = dataImport(AccData)
## names(segData)
}
