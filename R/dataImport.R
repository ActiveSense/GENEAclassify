

#' Loads the data into R and creates format required for segmentation.
#' 
#' @title Data import function
#' @param binfile File path to binary data to be segmented.
#' @param downsample Rate to downsample the data, defaults to every 100th observation. For no downsampling set NULL.
#' @param start Where to start reading observations.
#' @param end Where to end reading observations.
#' @param Use.Timestamps To use timestamps as the startand end time values this has to be set to TRUE. (Default FALSE)
#' @param radians calculate degrees rotation in radians.
#' @param keep_raw_data Keep the raw data for calculating steps using stepcounter. 
#' @param mmap.load Default is (.Machine$sizeof.pointer >= 8). see mmap for more details
#' @param ... additional arguments passed through.
#' @details Reads in the binary data file and extracts the information required for the segmentation procedure.
#' @return Returns a list containing a matrix of the data including the x, y and z axis data, vectors of the up down (elevation) 
#' and degrees (rotation), a vector of time stamps, a vector of vector magnitudes and the serial number of the device. 
#' @export
#' @import GENEAread
#' @examples
#' ##   segData = dataImport(bindata = file.path(system.file(package = "GENEAread"),
#' ##                                                         "binfile",
#' ##                                                         "TESTfile.bin"))
#' ## 
#' ## segData1 = dataImport(AccData)
#' ## names(segData)

dataImport = function(binfile, 
                      downsample = 100, 
                      start = NULL, 
                      end = NULL, 
                      Use.Timestamps = FALSE, 
                      radians = FALSE,
                      keep_raw_data = TRUE,
                      mmap.load = (.Machine$sizeof.pointer >= 8), pagerefs, ...) {

  # Note to bring everything to binfile name change. - To be removed
  if (missing(binfile)){stop("bindata has been renamed as binfile. Please rename variable")}

  # This should be an option!
  if (keep_raw_data){
    binaryDataOut = read.bin(binfile = binfile,
                             start = start,
                             end = end, 
                             Use.Timestamps = Use.Timestamps, 
                             mmap.load = mmap.load,
                             calibrate = TRUE,
                             downsample = 1,
                             pagerefs = pagerefs)
    
    RawData = binaryDataOut$data.out
  }
  
  if (is.null(downsample)) {
      
    binaryData = binaryDataOut
      
    binaryDataOut = binaryData$data.out
      
  } else {
    binaryData = read.bin(binfile, 
                          start = start, 
                          end = end, 
                          Use.Timestamps = Use.Timestamps, 
                          mmap.load = mmap.load,
                          calibrate = TRUE, 
                          downsample = downsample,
                          pagerefs = pagerefs)
    
    binaryDataOut = binaryData$data.out
   }
   
  if (keep_raw_data){
    binaryDataOut = RawData
  }
   
  serial = binaryData$header["Device_Unique_Serial_Code", ][[1]]

  rightWrist = grepl("right wrist", binaryData$header["Device_Location_Code", ][[1]])
  leftWrist = grepl("left wrist", 
            binaryData$header["Device_Location_Code", ][[1]]) | grepl("[[:blank:]]", 
            gsub("", " ", binaryData$header["Device_Location_Code", ][[1]]))

  if(!leftWrist & !rightWrist){
      warning("Note: data assumed to be left wrist")
  }

  Intervals = get.intervals(binaryData, 
                            length = NULL, 
                            incl.date = TRUE,
                            size = 1e6)

  ## Extract time, light and temp data
  Time = Intervals[, "timestamp"]
  Light = binaryData$data.out[, "light"]
  Temp = binaryData$data.out[, "temperature"]

  rm(binaryData)

  ## extract the up/down and rotation data
  dataUpDown = updown(Intervals)
  dataDegrees = degrees(Intervals)

  if (rightWrist) {
      dataUpDown = -1 * dataUpDown
  }

  vecMagnitude = abs(sqrt(rowSums((Intervals[, c("x", "y", "z")])^2)) - 1)
    
  if (radians){
      
      dataRadians = radians(Intervals)
      
      geneaBin = list(Data = Intervals, 
                      Freq = downsample, 
                      UpDown = dataUpDown, 
                      Degrees = dataDegrees, 
                      Radians = dataRadians,
                      Time = Time, 
                      Light = Light, 
                      Temp = Temp, 
                      Magnitude = vecMagnitude, 
                      RawData = binaryDataOut, 
                      Serial = serial)
    } else {
      geneaBin = list(Data = Intervals, 
                      Freq = downsample, 
                      UpDown = dataUpDown, 
                      Degrees = dataDegrees, 
                      Time = Time, 
                      Light = Light, 
                      Temp = Temp, 
                      Magnitude = vecMagnitude, 
                      RawData = binaryDataOut, 
                      Serial = serial)
    }
    
    class(geneaBin) = c(class(geneaBin), "GENEAbin")
    
    return(geneaBin)
}


#' Extract up/down time series
#' 
#' Extract code from \code{positionals} to perform data conversion to up/down time series.
#' Input is expected to be result of \code{get.intervals}.
#' @title Extract data relating to the up/down component
#' @param x data output from \code{get.intervals}
#' @return The up/down vertical elevation data (y-axis)
#' @export
#' @keywords internal
#' @examples
#'     d1 = matrix(c(100, 101, -0.79, -0.86, -0.17, -0.14, 0.53, 0.46), 
#'         nrow = 2, ncol = 4)
#'     colnames(d1) = c("timestamp", "x", "y", "z")
#'     updown(x = d1)

updown = function(x) {
    
    numerator = x[, "y"]
    
    # magnitude
    denominator = sqrt(rowSums(x[, c("x", "y", "z")]^2)) 
    
    ud = (-acos(numerator / denominator) * 180 / pi + 90)
    
    return(ud)
}

#' Extract data relating to the rotation component.
#' 
#' Called by \code{dataImport}.
#' Note: the "+ 1" has been removed from the original implementation.
#' @title Extract rotation time series
#' @param x data output from get.intervals
#' @return The degrees (rotation) data.
#' @export
#' @keywords internal
#' @examples
#'    d1 = matrix(c(100, 101, -0.79, -0.86, -0.17, -0.14, 0.53, 0.46), 
#'         nrow = 2, ncol = 4)
#'    colnames(d1) = c("timestamp", "x", "y", "z")
#'    degrees(x = d1)

degrees = function(x) {
    
    magnitude = sqrt(rowSums(x[, c("x", "z")]^2))
    
    deg = 361 * (sign(-x[, "x"]) * 180 * acos(-x[, "z"] / magnitude) / pi + 180) / 360
    
    deg = floor(deg) - 45 
    
    deg = ifelse(deg < 0, deg + 360, deg)
    
    return(deg)
}

#' Extract data relating to the rotation component.
#' 
#' Called by \code{dataImport}.
#' Note: the "+ 1" has been removed from the original implementation.
#' @title Extract rotation time series in radians 
#' @param x data output from get.intervals
#' @return The degrees (rotation) data.
#' @export
#' @keywords internal
#' @examples
#'    d1 = matrix(c(100, 101, -0.79, -0.86, -0.17, -0.14, 0.53, 0.46), 
#'         nrow = 2, ncol = 4)
#'    colnames(d1) = c("timestamp", "x", "y", "z")
#'    degrees(x = d1)

radians = function(x) {
  
  magnitude = sqrt(rowSums(x[, c("x", "z")]^2))
  
  rad = sign(-x[, "x"]) * acos(-x[, "z"] / magnitude) + pi  

  return(rad)
}

