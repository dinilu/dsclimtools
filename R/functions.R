#' Load a cropped version of dsclim (downscaled TraCE21ka) .nc files
#'
#' @param file A string with full name (path and file name) of the .nc file to be loaded from the hard disk drive.
#' @param sf A sf object from the sf package to be used as template or mask to crop the downscaled trace .nc file.
#' @param proxy A logical value indicating whether the data should be loaded as a stars_proxy object (TRUE) or as a stars object (FALSE).
#'
#' @return A stars object with the downscaled trace data.
#' @export
#'
#' @import stars sf
.read_dsclim <- function(file, sf = NULL, proxy) {
  data <- stars::read_stars(file, proxy = proxy)
  if(!is.null(sf)){
    # data <- sf::st_crop(data, sf)
    data <- stars::st_extract(data, sf)
  }
  data
}



#' Load dsclim (downscaled TraCE21ka) files by years.
#'
#' @param folder A character string with the path to the folder where the dsclim (downscaled TraCE21ka) files, ordered in folders by variables.
#' @param var A character string with the name of the variables to be loaded.
#' @param y_start A number with the first year to be loaded. Years are provided as cal BP (1950 calendar year = 1 cal BP)
#' @param y_end A number with the last year to be loaded.
#' @param rcp When loading future data (1991-2100), the Representative Concentration Pathway that want to be loaded. Valid values are: "rcp2.6", "rcp4.5", "rcp6.0", or "rcp8.5".
#' @param gcm When loading future data (1991-2100), the General Circulation Models that want to be loaded. Valid values are: "CESM1-CAM5", "CSIRO-Mk3-6-0", or "IPSL-CM5A-MR".
#' @param calendar_dates A logical value indicating whether dates should be corrected to calendar dates in the output.
#' @param sf A sf object to be used for cropping the object, can be points or polygons.
#' @param proxy A logical value indicating whether the data should be loaded as a stars_proxy object (TRUE) or as a stars object (FALSE).
#'
#' @return A stars object with the dsclim (downscaled TraCE21ka) data, cropped using the sf object if provided.
#'
#' @details Years should be in the format calibrated Before Present (e.g. 0 for loading year 1950 in the gregorian calendar).
#' @export
#'
#' @examples
#' \dontrun{
#' read_dsclim("data/dsclim/",
#'             var = "tasmax",
#'             -22000,
#'             -19950,
#'             proxy = FALSE)
#' read_dsclim("data/dsclim/", var = "tasmax",
#'             41,
#'             150,
#'             rcp = "rcp2.6",
#'             gcm = "CESM1-CAM5",
#'             calendar_dates = TRUE,
#'             proxy = FALSE)
#' }
read_dsclim <- function(folder, var, y_start, y_end, rcp = NULL, gcm = NULL, calendar_dates = FALSE, sf = NULL, proxy = TRUE){

  if(proxy & calendar_dates){
    stop("Dates shouldn't be changed in a stars_proxy object to avoid errors when loading
         the data from files in the deferred read/write loop.")
  }

  y_start <- as.numeric(y_start)
  y_end <- as.numeric(y_end)

  if(y_start == 0){
    y_start <- -1
  }

  if(y_end == 0){
    y_end <- -1
  }

  if(y_start >= y_end){
    stop("Starting year same or higher than ending year. Provide a y_start value lower than y_end.")
  }

  if(y_end > 40){
    if(is.null(rcp) || is.null(gcm)){
      stop("To load future data a Representative Concentration Pathway (rcp) and a Global Circulation Model (gcm) has to be specified.")
    }
  }

  # if(y_start < 0 & y_end > 0){
  #   y_seq <- c(y_start:-1, 1:y_end)
  # } else {
  #   y_seq <- c(y_start:y_end)
  # }
  y_seq <- y_start:y_end
  y_seq <- y_seq[which(y_seq != 0)]
  y_seq <- split(y_seq, y_seq > 40)

  # files <- paste0(folder, "/TraCE21ka/", var, "/", var, y_seq, ".nc")
  files <- NULL
  if(!is.null(y_seq$'FALSE') & !is.null(y_seq$'TRUE')){
    stop("You requested to load data that include past and future data.
         Please specify a time period that include past or future data only.")
  }
  if(!is.null(y_seq$'FALSE')){
    files <- append(files, paste0(folder, "/trace21ka/", var, "/", var, y_seq$'FALSE', ".nc"))
  }
  if(!is.null(y_seq$'TRUE')){
    files <- append(files, paste0(folder, "/", rcp, "/", gcm, "/", var, "/", var, y_seq$'TRUE', ".nc"))
  }

  if(length(files) > 1){
    data <- lapply(files, FUN=.read_dsclim, sf = sf, proxy = proxy)
    data <- do.call(c, data)
  } else {
    data <- .read_dsclim(files, sf = sf, proxy = proxy)
  }

  names(data) <- var

  if(calendar_dates){
    if(sign(y_start) == 1 && sign(y_end) == 1){
      data <- time_2_calendar_dates(data, y_start, y_end)
    }
    if(sign(y_start) == -1 && sign(y_end) == 1) {
      data <- time_2_calendar_dates(data, y_start+1, y_end)
    }
    if(sign(y_start) == -1 && sign(y_end) == -1){
      data <- time_2_calendar_dates(data, y_start+1, y_end+1)
    }
  }

  data
}



#' Create a calendar dates by specified time intervals.
#'
#' @param y_start A number with the first year for the required dates.
#' @param y_end A number with the last year for the required dates.
#' @param by A number or string (as in lubridate package) to be used as interval to get the calendar dates. Default to "1 month" which is returning monthly dates (first day of the month).
#'
#' @return A vector with dates from the y_start to the y_end years at the time intervals specified by "by" argument.
#' @export
#'
#' @examples
#' \dontrun{
#' calendar_dates(41, 150, by = "1 month")
#' calendar_dates(41, 150, by = "1 year")
#' calendar_dates(-22000, 150, by = "1 year")
#' }
calendar_dates <- function(y_start, y_end, by = "1 month"){
  cal_start <- lubridate::ymd("1950-01-01") + lubridate::years(y_start)
  cal_end <- lubridate::ymd("1950-12-01") + lubridate::years(y_end)
  seq(cal_start, cal_end, by = by)
}



#' Modify time to calendar years
#'
#' Modify a stars object to change the time dimension to calendar years between specified starting and ending years.
#'
#' @param data A stars object with the data to be modified.
#' @param y_start A number with the starting year of the data (in calibrated Before Present format).
#' @param y_end A number with the ending year of the data (in calibrated Before Present format).
#' @param by A number or string specifying the interval of the dates to be used in the new stars object.
#'
#' @return A stars object as in data argument but with changed time dimension.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read_dsclim("data/dsclim/",
#'                     var = "tasmax",
#'                     41,
#'                     150,
#'                     rcp = "rcp2.6",
#'                     gcm = "CESM1-CAM5",
#'                     calendar_dates = FALSE,
#'                     proxy = FALSE)
#' time_2_calendar_dates(data,
#'                       41,
#'                       150,
#'                       by = "1 month")
#' }
time_2_calendar_dates <- function(data, y_start, y_end, by = "1 month") {
  cal_dates <- calendar_dates(y_start, y_end, by = by)
  # stars::st_dimensions(data)$time$values <- cal_dates
  data <- stars::st_set_dimensions(data, which = "time", values = cal_dates)
  data
}




#' Read original TraCE21ka datasets
#'
#' Read original TraCE21ka data, returning a stars object. This function differs from the ones in the dsclim package, which return the data in the format of the climate4r framework.
#'
#' @param folder character string indicating the path to the folder with the TraCE21ka dataset.
#' @param var character string indicating the name of the variable that is to be readed.
#' @param sf a sf vector object for which the TraCE21ka data should be limited (cropped).
#'
#' @return The function returns an stars object, usually in the form of raster data, except if the argument sf is provided, in which case the object is a vector stars object.
#' @export
#'
#' @examples
#' \dontrun{
#' read_trace("data/TraCE21ka/", "TSMX")
#' }
read_trace <- function(folder, var, sf = NULL){
  # folder <- "../../Public/Data/TraCE21ka/"
  # var <- "PRECC"
  # sf <- st_set_crs(points[1,], 4326)
  files <- list.files(paste0(folder, "/", var), full.names = TRUE, pattern=".nc")
  paste(files)
  data <- lapply(files, FUN = stars::read_ncdf, proxy = FALSE)
  # data <- lapply(files, FUN = stars::read_stars) # Doesn't work because uneven grid.
  data <- do.call(c, data)

  data <- time_2_calendar_dates(data, -21999, 40)

  if(var %in% c("TS", "TSMX", "TSMN")){
    data <- kelvin2celsius(data)
  }
  if(var == "PRECC"){
    data <- flux2mm(data)
  }

if(!is.null(sf)){
    lat <- st_coordinates(sf)[,2]
    lats <- stars::st_dimensions(data)$lat$values
    i <- findInterval(lat, lats)

    data <- dplyr::filter(data, lat > floor(lats[i]), lat < ceiling(lats[i+1]))
    # data <- sf::st_crop(data, sf)
    data <- stars::st_extract(data, at = sf)
  }
  data
}


#' Degrees Kelvin to Celsius
#'
#' Transform stars objects from degrees Kelvin to degrees Celsius
#'
#' @param data A stars object with temperature in degrees Kelvin.
#'
#' @return A stars object with the temperature in degrees Celsius
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read_trace("data/TraCE21ka/", "TSMX")
#' kelvin2celsius(data)
#' }
kelvin2celsius <- function(data) {
  var <- names(data)
  data_class <- class(data[[var]])
  data[[var]]  <- as.numeric(data[[var]]) - 273.15
  class(data[[var]]) <- data_class
  units <- list(numerator = "\u00B0C", denominator = character(0))
  class(units) <- "symbolic_units"
  attr(data[[var]], "units") <- units
  data
}


#' Flux 2 millimetters
#'
#' Transform stars objects from precipitation in flux to precipitation in millimeters
#'
#' @param data A stars object with precipitation in flux
#'
#' @return A stars object with precipitation in millimeters
#' @export
#'
#' @examples
#' \dontrun{
#' data <- read_trace("data/TraCE21ka/", "PRECC")
#' flux2mm(data)
#' }
flux2mm <- function(data){
  var <- names(data)
  data_class <- class(data[[var]])
  data[[var]]  <- as.numeric(data[[var]]) * 2592000000
  class(data[[var]]) <-data_class
  units <- list(numerator = "mm", denominator = character(0))
  class(units) <- "symbolic_units"
  attr(data[[var]], "units") <- units
  data
}

