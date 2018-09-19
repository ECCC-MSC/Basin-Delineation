#===============================================================================
#' @title Create spatial object from Hydat connection
#'
#' @description  Create a SpatialPointsDataFrame from the "STATIONS" Hydat table.
#'
#' @param con  An open SQLite database connection to the HYDAT database or a
#' character string pointing to
#' a SQLite database
#'
#' @param ...  Optional character vector of stations to subset (passed to
#' HYDAT::StationMetaData))
#'
#' @return A SpatialPointsDataFrame for the Hydat stations
#'
#' @export
#===============================================================================
SpatialHydat <- function(con, ...){
  if (class(con) == "character"){
    if (!requireNamespace("RSQLite")){
      print("Missing library: 'RSQLite'. Please install library")
      return()
    }
    con <- RSQLite::dbConnect(RSQLite::SQLite(), con)
    }
  if (!requireNamespace("HYDAT")){
    print("Missing library: 'HYDAT'. Please install library")
    return()
    }
  Hyd <- HYDAT::StationMetadata(con, ...)
  coord <- Hyd[, c("longitude", "latitude")]
  output <- sp::SpatialPointsDataFrame(coords = coord, data = Hyd,
        proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
  return(output)
}

#===============================================================================
#' @title Spatial ECDE
#'
#' @description  Create a SpatialPointsDataFrame from ECDE table export.
#'
#' @param data.file  A file exported from ECDE
#'
#' @param data.format  (optional) Character string describing file type.
#' Options include: ('csv', 'tb0', 'dbf' or 'pt2').  If missing, the file
#' extension of data.file will be used as a guess.
#'
#' @return A SpatialPointsDataFrame for the Hydat stations
#'
#' @export
#===============================================================================
SpatialECDE <- function(data.file, data.format){
  # Read in data and deal with any formatting issues
  if (missing(data.format)){
    data.format <- gsub("^.*\\.(.+)$", "\\1", data.file)
  }

  if (data.format == "csv") {
    input <- read.csv(data.file, stringsAsFactors = F)
      if ("DISCLAIMER" %in% input[, 1]){
        headers <- names(input)
        input <- read.csv(data.file, skip = 2, header = F, blank.lines.skip = T,
                          stringsAsFactors = F)

        names(input) <- headers
        input <- input[-c(nrow(input), nrow(input) - 1), ] # remove DISCLAIMER rows
      }
    # clean leading spaces
    i <- lapply(input, class) == "character"
    input[, i] <- sapply(input[, i], gsub, pattern = "^ ", replacement = "")
  }else if (data.format == "dbf"){
    if (!requireNamespace("foreign")){
      print("Missing library: 'foreign', try another format or install package")
    }

    input <- foreign::read.dbf(data.file)
  }else if (data.format == "tb0"){
    input <- read.tb0(data.file)
  }else if (data.format == "pt2") {
    input <- read.pt2(data.file)
  }else {
    print("Could not read file. Try specifying a data format")
  }

  if (data.format != "dbf") {
    input[] <- lapply(input, bool_check)
  }

  input[, sapply(input, class) == "factor"] <- lapply(input[, sapply(input, class) == "factor"], as.character)

  # Change names so that they match HYDAT.  (specifically what is returned by hydat R package)
  names(input) <- tolower(names(input))  # Hydat R package has lowercase convention
  names(input) <- sub("^station$",      "station_number",      names(input))
  names(input) <- sub("^stationname$",  "station_name",        names(input))
  names(input) <- sub("^hydstatus$",    "hyd_status",          names(input))
  names(input) <- sub("^prov$",         "prov_terr_state_loc", names(input))
  names(input) <- sub("^drainagearea$", "drainage_area_gross", names(input))
  names(input) <- sub("^prov$",         "prov_terr_state_loc", names(input))
  names(input) <- sub("^region$",       "regional_office",     names(input))
  names(input) <- sub("^sed$",          "sed_status",          names(input))
  names(input) <- sub("^realtime$",     "real_time",           names(input))

  coord <- input[, c("longitude", "latitude")]
  output <- sp::SpatialPointsDataFrame(coords = coord, data = input,
                proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )

  return(output)
}

#==============================================================================
#' @title Create SpatialPointsDataFrame from csv
#'
#' @description Creates an R SpatialPointsDataFrame from a csv
#'
#' @param file path to *.csv file
#'
#' @param ID character, column name giving a unique ID for each entry. If Hydat 'station_number' entries are
#' available for every station, this is the recommended column to use.
#'
#' @param header_x character, column name containing x coordinate information
#'
#' @param header_y character, column name containing y coordinate information
#'
#' @param CRS (optional) character, proj4 string specifying the projection
#' information of the points. If missing, the default is WGS84
#'
#' @return SpatialPointsDataFrame
#'
#' @export
#==============================================================================
SpatialCSV <- function(file, header_x, header_y, CRS, ID){
  if (missing(CRS)){
    CRS <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  }

  CRS <- sp::CRS(CRS)
  file_data <- read.csv(file, stringsAsFactors = F)

  if (!missing(ID)){
    if ("station_number" %in% names(file_data) &  "ID" != "station_number"){
      warning("station_number column has been overwritten")
    }
    file_data$station_number <- file_data[, ID]
  }else{
    warning("Missing 'station_number' column. Unexpected behaviour may result with basin delineation tools")
  }

  coords <- file_data[, c(header_x, header_y)]
  output <- sp::SpatialPointsDataFrame(coords, data = file_data, proj4string = CRS)

  return(output)
}

#===============================================================================
#' @title Read pt2 file
#'
#' @description  Read in a pt2 file from EC Data explorer output table. Does
#' some type conversion to make the output the same as the read.tb0 and read.dbf
#' functions for ECDE
#'
#' @param x  pt2 file from EC Data explorer export
#'
#' @return A data.frame with (hopefully) appropriate data types
#'
#' @export
#===============================================================================
read.pt2 <- function(x){
  # read in header information
  n.attributes <- read.table(x, skip = 15, nrows = 1)[1, 2]
  names_attributes <- readLines(x, n = n.attributes * 2 + 16)[17:(n.attributes * 2 + 16)]
  data.types <- names_attributes[seq(2, length(names_attributes), 2)]

  names_attributes <- names_attributes[seq(1, length(names_attributes), 2)]
  names_attributes <- matrix(unlist(strsplit(names_attributes, " ")), ncol = 3, byrow = T)[, 3]

  # read in data
  out <- read.table(x, skip = (n.attributes * 2) + 16 + 2, stringsAsFactors = F)
  out <- out[, -c(1, 2)]
  names(out) <- names_attributes

  # Some things are stored as factors and they get screwed up during read-in. Correct this:

  # format text strings
  which_factors <- grepl(pattern = ".* oneof", data.types)
  factor_defs <- gsub("^.* oneof ", "",  data.types[which_factors])
  factor_defs <- gsub("\" \"",      ",", factor_defs)
  factor_defs <- gsub("[\"]",       "",  factor_defs)
  factor_defs <- gsub("^ ",         "",  factor_defs)

  # break up text strings
  factor_defs <- strsplit(factor_defs, ",")
  factor_defs <- lapply(factor_defs, strsplit, split = "=")
  names(factor_defs) <- names(out)[which_factors]

  # replace digit indices with factors for each factor variable, then convert to character
  for (i in seq(1, length(data.types))){
    if (which_factors[i] == TRUE){
      df <- matrix( (unlist(factor_defs[names(out)[i]])), ncol = 2, byrow = T)
      replacement <- factor(out[, i], levels = df[, 2], labels = df[, 1])
      out[, i] <- replacement
    }
  }
  return(out)
}

#===============================================================================
#' @title Read *.tb0
#'
#' @description  Read in a tb0 file from EC Data explorer output table. Does some type
#' conversion to make the output the same as the read.pt2 and read.dbf functions for ECDE
#'
#' @param x  tb0 file from EC Data explorer export
#'
#' @return A data.frame with (hopefully) appropriate data types
#'
#' @export
#===============================================================================
read.tb0 <- function(x){
  colnam <- read.table(x, skip = 26, nrows = 1, stringsAsFactors = F)
  colnam <- colnam[1, -1]
  data.cols <- read.table(x, skip = 33, stringsAsFactors = F)
  names(data.cols) <- colnam
  out <- data.cols
  return(out)
}

#===============================================================================
#' @title Convert vector to logical
#'
#' @description  Changes vectors into Boolean when appropriate
#'
#' @param vec  A vector (any type)
#'
#' @return the original unchanged vector, or if possible a boolean
#'
#' @keywords internal
#===============================================================================
bool_check <- function(vec){
  if (class(vec) == "character") {
    vec.uniq <- trimws(tolower(unique(vec)))
    if ( ("true" %in% vec.uniq | "false" %in% vec.uniq) & (length(vec.uniq) <= 2)){
      vec <- as.logical(trimws(vec))
    }
  }
  return(vec)
}
