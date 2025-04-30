
#' Function to calculate the number and variance of the steps in the data.
#'
#' @title Step Counter
#' @param AccData The data to use for calculating the steps. This should either an AccData object or a vector.
#' @param samplefreq The sampling frequency of the data, in hertz,
#' when calculating the step number (default 100).
#' @param filterorder single integer, order of the Chebyshev bandpass filter,
#' passed to argument n of \code{\link[signal]{cheby1}}.
#' @param boundaries length 2 numeric vector specifying lower and upper bounds
#' of Chebychev filter (default \code{c(0.5, 5)} Hz),
#' passed to argument W of \code{\link[signal]{butter}} or \code{\link[signal]{cheby1}}.
#' @param Rp the decibel level that the cheby filter takes, see \code{\link[signal]{cheby1}}.
#' @param plot.it single logical create plot of data and zero crossing points (default \code{FALSE}).
#' @param hysteresis The hysteresis applied after zero crossing. (default 100mg)
#' @param fun character vector naming functions by which to summarize steps.
#' "count" is an internally implemented summarizing function that returns step count.
#' @param verbose single logical should additional progress reporting be printed at the console? (default FALSE).
#' @return Returns a vector with length fun.
#' @export
#' @importFrom signal cheby1
#' @examples
#' d1 <- sin(seq(0.1, 100, 0.1))/2 + rnorm(1000)/10 + 1
#' Steps4 = stepCounter(d1)
#' length(Steps4)
#' mean(Steps4)
#' sd(Steps4)
#' plot(Steps4)

.timing_env <- new.env()

.timing_env$start_total <- Sys.time()
.timing_env$last_step <- Sys.time()

trace_time <- function(message) {
  current_time <- Sys.time()
  
  # Calculate time since the last trace_time call (or script start for the first call)
  step_duration <- current_time - .timing_env$last_step
  
  # Calculate total time since the initial setup
  total_duration <- current_time - .timing_env$start_total
  
  # Print the timing information neatly
  cat(sprintf("--- Timing Trace --- %s\n", message))
  # Use as.numeric to get the duration value, attr for the units (usually "secs")
  cat(sprintf("  Step Duration: %.3f %s\n", as.numeric(step_duration), attr(step_duration, "units")))
  cat(sprintf("  Total Elapsed: %.3f %s\n", as.numeric(total_duration), attr(total_duration, "units")))
  cat("--------------------\n")
  
  # Update the last step time for the next calculation
  .timing_env$last_step <- current_time
}

stepCounter <- function(AccData, 
                        samplefreq = 100,
                        filterorder = 2,
                        boundaries = c(0.5, 5), # 
                        Rp = 3,
                        plot.it = FALSE, 
                        hysteresis = 0.05, 
                        verbose = verbose,
                        fun = c("GENEAcount", "mean", "sd", "mad")) {
  
  if (missing(AccData)) { stop("data is missing") }
  if (!is.character(fun)) { stop("fun must be character vector of function names") }
  if (length(fun) < 1L) { stop("fun must name at least one function") }
  
  # Going to remove the amplitude variable from this step counter. Can come back to this later. 
  res <- numeric(length(fun))
  names(res) <- fun
  
  #### Check whether an AccData object or a vector ####
  
  # If Not an AccData object is it a numerical vector (Can be timestamps and a vector)
  if (inherits(AccData, "AccData")){
    StepData = AccData$data.out[,3]
    samplefreq = AccData$freq
  } else if (inherits(AccData, "numeric")){
    StepData = AccData
    if (missing(samplefreq)){
      warning("No samplefreq is given. samplefreq set to default, samplefreq = 100")
      samplefreq = 100
    }
  } else if (inherits(AccData,"matrix") & dim(AccData)[2] == 2 | 
             inherits(AccData, "data.frame") & dim(AccData)[2] == 2){
    StepData = AccData[,2]
    if (missing(samplefreq)){
      warning("No samplefreq is given. samplefreq set to default, samplefreq = 100")
      samplefreq = 100
    }
  } else {
    stop("Step Counter must use either an AccData object, Numerical Vector or a 2D Matrix of time and StepData")
  }
  
  Filter <- cheby1(n = filterorder,                               # order of filter
                   Rp = Rp,                                       # ripple of band pass
                   W = boundaries/samplefreq,                     # lower then upper frequencies of band pass
                   type = "pass",
                   plane = "z")
  
  #### Apply the band pass filter ####
  filteredData = signal::filter(Filter, StepData) 
  
  # ==================================
  # NEW IMPLEMENTATION
  # ==================================
  
  samples <- length(filteredData)
  
  # 1. Determine the zone for each point: -1 (below), 0 (between), 1 (above)
  #    Using integers (1L, -1L, 0L) can be slightly more efficient
  zone <- ifelse(filteredData > hysteresis, 1L, 
                 ifelse(filteredData < -hysteresis, -1L, 0L))
  
  # 2. Find where the zone changes compared to the previous point
  #    Create a 'lagged' version of the zone vector. Assume starts in zone 0 or -1 
  #    (similar to original loop starting with state = -1, let's assume 0 is safe default)
  prev_zone <- c(0L, zone[-samples]) 
  
  # 3. Identify the exact indices where a transition *into* the target zone occurs
  #    Enters zone 1 (crossed +hysteresis upwards)
  enters_pos_zone_idx <- which(zone == 1L & prev_zone != 1L)
  #    Enters zone -1 (crossed -hysteresis downwards)
  enters_neg_zone_idx <- which(zone == -1L & prev_zone != -1L)
  
  # 4. Combine and sort all relevant transition events by their index
  if (length(enters_pos_zone_idx) > 0 || length(enters_neg_zone_idx) > 0) {
    
    events <- data.frame(
      idx = c(enters_pos_zone_idx, enters_neg_zone_idx),
      type = rep(c(1L, -1L), c(length(enters_pos_zone_idx), length(enters_neg_zone_idx)))
    )
    events <- events[order(events$idx), , drop = FALSE]
    
    # 5. Find the 'up' events (type 1) that immediately follow a 'down' event (type -1)
    #    This implements the hysteresis logic: count 'up' only after a 'down'
    prev_event_type <- c(0L, events$type[-nrow(events)]) # Type of the *previous event*
    
    # Indices of the valid step start points (upward crossing after a downward one)
    valid_step_start_indices <- events$idx[events$type == 1L & prev_event_type == -1L]
    
    # 6. Calculate durations between consecutive valid step starts
    if (length(valid_step_start_indices) > 1) {
      cadence_samples <- diff(valid_step_start_indices)
    } else {
      cadence_samples <- numeric(0) # No complete steps found
    }
    
  } else {
    # No threshold crossings at all
    valid_step_start_indices <- numeric(0)
    cadence_samples <- numeric(0)
  }
  
  # 7. Convert durations to seconds
  cadence <- cadence_samples / samplefreq
  
  # ==================================
  # OLD IMPLEMENTATION (TOO SLOW)
  # ==================================
  
  # state = -1                                                       # initialise step state
  # interval = 0                                                     # initialise the interval counter
  # cadence = numeric(0)                                             # initialise first element of array for intervals
  # samples = length(filteredData)                                   # loop through all samples
  # 
  # for (a in 1:samples) {
  #   if ((filteredData[a] > hysteresis) && (state < 0)){            # new step started
  #     state = 1                                                    # set the state
  #     cadence[length(cadence)] = interval + 1                      # write the step interval
  #     cadence[length(cadence) + 1] = 0                             # initialise to record the next step    
  #     interval = 0                                                 # reset the step counter
  #   } else if ((-1*filteredData[a] > hysteresis) && (state > 0)) { # hysteresis reset condition met
  #     state = -1                                                   # reset the state
  #     interval = interval + 1                                      # increment the interval
  #   } else {
  #     interval = interval + 1                                      # increment the interval
  #   }
  #   cadence[length(cadence)] = interval                            # capture last part step
  #   
  #   print(paste0("--- INTERVAL --- # ", a))
  #   print(cadence)
  #   
  #   print(paste0("--- CADENCE --- # ", a))
  #   print(cadence)
  # }
  # 
  # cadence = cadence/samplefreq                                     # divide by the sample frequency to get seconds
  
  if ("GENEAcount" %in% fun) {
    res["GENEAcount"] <- length(cadence)
    fun <- fun[fun != "GENEAcount"]
  }
  
  if ("mean" %in% fun) {
    if (length(cadence) < 2){
      res["mean"] <- 0
      fun <- fun[fun != "mean"]
    } else {
      res["mean"] <- 60 / mean(cadence, na.rm = T)
      fun <- fun[fun != "mean"]
    }
  }
  
  if ("median" %in% fun) {
    if (length(cadence) < 2){
      res["median"] <- 0
      fun <- fun[fun != "median"]
    } else {
      res["median"] <- 60 / median(cadence, na.rm = T)
      fun <- fun[fun != "median"]
    }
  }
  
  for (i in fun) {
    val <- try(get(x = i, mode = "function")(diff(cadence)))
    if (is(val, class2 = "try-error")) { val <- 0 }
    if (is.na(val)) { val <- 0 }
    res[i] <- val
  }
  
  return(res)
  
}

#' Peak detect function 
#' 
#' @title Find_Peaks
#' @param x timeseries signal 
#' @param m The number of points either side of the peak to required to be a peak. 
#' @details Finds the peaks and valleys within the signal passed to the function.  
#' @keywords Internal
#' 

find_peaks <- function(x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  return(pks)
}



