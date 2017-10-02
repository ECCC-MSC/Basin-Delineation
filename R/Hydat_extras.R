#'  SpatialHydat
#'
#' @title Create a SpatialPointsDataFrame from a Hydat database connection.
#' @description  Create a SpatialPointsDataFrame from the "STATIONS" Hydat table.
#' @param con  An open SQLite database connection to the HYDAT database
#' @param ...  Optional character vector of stations to subset (passed to HYDAT::StationMetaData))
#' @return A SpatialPointsDataFrame for the Hydat stations
#' @export
SpatialHydat <- function(con, ...){
  # con <- dbConnect(RSQLite::SQLite(), "path to hydat sqlite database") 
  Hyd <- HYDAT::StationMetadata(con, ...)
  coord <- Hyd[,c("longitude","latitude")]
  #Hyd <- Hyd[,-which(names(a) %in% c("longitude", "latitude"))]
  output <- sp::SpatialPointsDataFrame(coords = coord, data=Hyd, proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
  return(output)
}

#'  SpatialECDE
#'
#' @title Create a SpatialPointsDataFrame from ECDE table export.
#' @description  Create a SpatialPointsDataFrame from ECDE table export.
#' @param data.file  A file exported from ECDE
#' @param data.format  (optional) Character string describing file type.  Options include: ('csv', 'tb0', 'dbf' or 'pt2').
#'   If missing, will use file extension to guess.
#' @return A SpatialPointsDataFrame for the Hydat stations
#' @export
SpatialECDE <- function(data.file, data.format){
  # Read in data and deal with any formatting issues
  if (missing(data.format)){
    data.format = gsub("^.*\\.(.+)$","\\1",data.file)
  }
  if (data.format == 'csv') {
    input <- read.csv(data.file, stringsAsFactors = F)
      if ("DISCLAIMER" %in% input[,1]){
        headers <- names(input)
        input <- read.csv(data.file, skip=2, header = F, blank.lines.skip = T, stringsAsFactors = F)
        names(input) <- headers
        input <- input[-c(nrow(input),nrow(input)-1),] # Get rid of DISCLAIMER rows
      }
    # clean leading spaces
    i <- lapply(input, class) == "character"
    input[,i] <- sapply(input[,i], gsub, pattern="^ ", replacement="")
  }else if (data.format == 'dbf'){
    if (!requireNamespace("foreign")){print("Missing library: 'foreign', try another format or install package")}
    input <- foreign::read.dbf(data.file)
  }else if (data.format == 'tb0'){
    input <- read.tb0(data.file)
  } else if (data.format == 'pt2') {
    input <- read.pt2(data.file)
  } else {
    print("Could not read file. Try specifying a data format")
  }
  if (data.format != "dbf") {
    input[] <- lapply(input, bool.check)
  }
  input[,sapply(input, class)=="factor"] <- lapply(input[,sapply(input, class)=="factor"], as.character)
  
  # Change names so that they match HYDAT.  (specifically what is returned by hydat R package)
  names(input) <- tolower(names(input))  # Hydat R package has lowercase convention
  names(input) <- sub("^station$", "station_number", names(input))
  names(input) <- sub("^stationname$", "station_name", names(input))  
  names(input) <- sub("^hydstatus$", "hyd_status", names(input)) 
  names(input) <- sub("^prov$", "prov_terr_state_loc", names(input)) 
  names(input) <- sub("^drainagearea$", "drainage_area_gross", names(input))  
  names(input) <- sub("^prov$", "prov_terr_state_loc", names(input)) 
  names(input) <- sub("^region$", "regional_office", names(input)) 
  names(input) <- sub("^sed$", "sed_status", names(input)) 
  names(input) <- sub("^realtime$", "real_time", names(input))
  
  coord <- input[,c("longitude","latitude")]
  output <- sp::SpatialPointsDataFrame(coords = coord, data=input, proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )
  return(output)
}

