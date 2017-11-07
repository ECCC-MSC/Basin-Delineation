#'  SpatialHydat
#'
#' @title Create a SpatialPointsDataFrame from a Hydat database connection.
#' @description  Create a SpatialPointsDataFrame from the "STATIONS" Hydat table.
#' @param con  An open SQLite database connection to the HYDAT database or a character string pointing to
#' a SQLite database
#' @param ...  Optional character vector of stations to subset (passed to HYDAT::StationMetaData))
#' @return A SpatialPointsDataFrame for the Hydat stations
#' @export
SpatialHydat <- function(con, ...){
  if (class(con)=="character"){
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
  coord <- Hyd[,c("longitude","latitude")]
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


#' Read *.pt2
#'
#' @title Read in a pt2 file from EC Data explorer output table
#' @description  Read in a pt2 file from EC Data explorer output table. Does some type
#' conversion to make the output the same as the read.tb0 and read.dbf functions for ECDE
#' @param x  pt2 file from EC Data explorer export
#' @return A data.frame with (hopefully) appropriate data types
#' @export
read.pt2 <- function(x){
  # read in header information
  n.attributes <- read.table(x, skip = 15, nrows = 1)[1,2]
  names.attributes <- readLines(x, n= n.attributes*2+16)[17:(n.attributes*2+16)]
  data.types <- names.attributes[seq(2,length(names.attributes),2)]
  names.attributes <- names.attributes[seq(1,length(names.attributes),2)]
  names.attributes <- matrix(unlist(strsplit(names.attributes,' ')), ncol=3, byrow = T)[,3]

  # read in data
  out <- read.table(x, skip=n.attributes*2+16+2, stringsAsFactors = F)
  out <- out[,-c(1,2)]
  names(out) <- names.attributes

  ## Note that some things are stored as factors and they get screwed up during read-in. Correct this:
  # this gets a little ugly...

  # format text strings
  which.factors <- grepl(pattern = '.* oneof', data.types)
  factor.defs <- gsub("^.* oneof ", '', data.types[which.factors])
  factor.defs <- gsub("\" \"", ',', factor.defs)
  factor.defs <- gsub("[\"]", "", factor.defs)
  factor.defs <- gsub("^ ", "", factor.defs)

  # break up text strings
  factor.defs <- strsplit(factor.defs, ',')
  factor.defs <- lapply(factor.defs, strsplit, split='=')
  names(factor.defs) <- names(out)[which.factors]

  # replace digit indices with factors for each factor variable, then convert to character
  for (i in seq(1,length(data.types))){
    if (which.factors[i] == TRUE){
      df <- matrix((unlist(factor.defs[names(out)[i]])), ncol=2, byrow = T)
      replacement = factor(out[,i], levels = df[,2], labels = df[,1])
      out[,i] <- replacement
    }
  }
  return(out)
}

#' Read *.tb0
#'
#' @title Read in a tb0 file from EC Data explorer output table
#' @description  Read in a tb0 file from EC Data explorer output table. Does some type
#' conversion to make the output the same as the read.pt2 and read.dbf functions for ECDE
#' @param x  tb0 file from EC Data explorer export
#' @return A data.frame with (hopefully) appropriate data types
#' @export
read.tb0 <- function(x){
  colnam <- read.table(x, skip = 26, nrows = 1, stringsAsFactors = F)
  colnam <- colnam[1,-1]
  data.cols <- read.table(x, skip = 33, stringsAsFactors = F)
  names(data.cols) <- colnam
  out <- data.cols
  return(out)
}

#'  Convert to Boolean
#'
#' @title Create a SpatialPointsDataFrame from ECDE table export.
#' @description  Changes vectors into Boolean when appropriate
#' @param vec  A vector (any type)
#' @return the original unchanged vector, or if possible a boolean
#' @keywords internal
#' @export
bool.check <- function(vec){
  if (class(vec) == "character") {
    vec.uniq <- trimws(tolower(unique(vec)))
    if (("true" %in% vec.uniq | "false" %in% vec.uniq) & (length(vec.uniq) <=2)){
      vec <- as.logical(trimws(vec))
    }
  }
  vec
}
