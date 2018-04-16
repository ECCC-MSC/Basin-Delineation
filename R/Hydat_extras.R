#===============================================================================
#' @title Create spatial object from Hydat connection
#'
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
    data.format = gsub("^.*\\.(.+)$","\\1",data.file)
  }
  if (data.format == 'csv') {
    input <- read.csv(data.file, stringsAsFactors = F)
      if ("DISCLAIMER" %in% input[,1]){
        headers <- names(input)
        input <- read.csv(data.file, skip=2, header = F, blank.lines.skip = T,
                          stringsAsFactors = F)

        names(input) <- headers
        input <- input[-c(nrow(input),nrow(input)-1),] # remove DISCLAIMER rows
      }
    # clean leading spaces
    i <- lapply(input, class) == "character"
    input[,i] <- sapply(input[,i], gsub, pattern="^ ", replacement="")
  }else if (data.format == 'dbf'){
    if (!requireNamespace("foreign")){
      print("Missing library: 'foreign', try another format or install package")
    }

    input <- foreign::read.dbf(data.file)
  }else if (data.format == 'tb0'){
    input <- read.tb0(data.file)
  } else if (data.format == 'pt2') {
    input <- read.pt2(data.file)
  } else {
    print("Could not read file. Try specifying a data format")
  }
  if (data.format != "dbf") {
    input[] <- lapply(input, BoolCheck)
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
  output <- sp::SpatialPointsDataFrame(coords = coord, data=input,
                proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )

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
  colnam <- colnam[1,-1]
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
BoolCheck <- function(vec){
  if (class(vec) == "character") {
    vec.uniq <- trimws(tolower(unique(vec)))
    if (("true" %in% vec.uniq | "false" %in% vec.uniq) & (length(vec.uniq) <=2)){
      vec <- as.logical(trimws(vec))
    }
  }
  return(vec)
}

#===============================================================================
#' @title Waterbody name from HYDAT station name
#'
#' @param station_name character string; station name from Hydat station
#'
#' @return character string corresponding to the hydrographic feature associated
#' with the station (i.e. the name of a lake, river, stream etc)
#'
#' @description Given a NHS / Hydat hydrometric station name, attempts to isolate
#' the name of the waterbody on which the station resides.
#'
#' @details The returned value can be used to query the canadian geographic
#'  names database.
#'
#' @return a character string of the waterbody name
#'
#' @examples
#' ParseStationName("PEARSON CREEK NEAR PROCTER")
#' ParseStationName("SCHIAVON CREEK BELOW HIGHEST DIVERSION NEAR THRUMS")
#' ParseStationName("MADAWASKA (RIVIER) EN AVAL DU BARRAGE TEMISCOUATA")
#'
#' @export
#===============================================================================
ParseStationName <- function(station_name){
  name <- toupper(station_name)

  join_terms_en <- c("AT", "ABOVE", "BELOW", "NEAR",
                     "EAST OF", "IN", "OUTFLOW", "BETWEEN", "WEST OF", "NORTH OF",
                     "SOUTH OF")
  join_terms_fr_1 <- c("EN AVAL", "EN AMONT", "CENTRALE", "PRES" )
  join_terms_fr_2 <- c("DU", "DE LA", "DE", "DE L'")
  join_terms_fr_3 <- c("AU", "A LA")

  #try to split on english join term
  splitters <- sapply(join_terms_en, pad, pads=' ')

  # find all splitting terms
  matches <- sapply(splitters,
                    function(p) gregexpr(p, text=station_name)[[1]][1])

  # if any there are any matches...
  if (!all(matches == -1)){
    matches[matches==-1] <- NA
    first_split <- splitters[which.min(matches)] # find first instance of splitting term
    name <- strsplit(station_name, split=first_split)[[1]][1]

  # if that fails, assume its a french name and move parenthetical description
  }else{
    name <- gsub(
      "([[:alnum:][:blank:]\\.\\'-]*) \\(([[:alnum:][:blank:]\\']*)\\).*",
      '\\2 \\1',
      station_name)
  }

  # Clean up name
  # remove anything after{ a comma, SITE, NO.
  name <- gsub("\\.", "", name) # strip off any periods
  name <- gsub("' ", "'",name) # get rid of spaces after apostrophes

  # deal with any remaining parenthetical material?

  return(name)
}

#===============================================================================
#' @title Geometry of waterbody from station name
#'
#' @desciption Attempts to figure out the geometry of a waterbody (point, polygon,
#' line) based on its name
#'
#' @param name name of waterbody
#'
#' @keywords internal
#'
#' @export
#===============================================================================
StationNameGeometry <- function(name){

  # List of common terms
  features_RIV_fr <- c("RIVIERE", "RIVIER", "CHENAL", "BRAS", "RUISSEAU",
                       "FLEUVE", "CANAL")
  features_LAC_fr <- c( "LAC", "BASSIN","ETANG", "BAIE", "LACS")
  features_LAC_en <- c( "LAKE", "BASIN","RESERVOIR", "POND", "SWAMP", "SLOUGH",
                        "SWAMP", "BAY", "POND", "LAKES" )
  features_RIV_en <- c("RIVER", "CREEK", "BROOK", "STREAM", "CHANNEL", "CANAL",
                       "RUN", "SPILLWAY", "TRIBUTARY", "DIVERSION", "FLOW")
  features_OTH_en <- c("DRAIN", "DITCH", "RILL", "SPRING", "SPRINGS", "COULEE")

  name <- toupper(name)
  name <- gsub("[^[:alnum:] ]", '', name) # remove wildcards if present

  # Try to determine expected shape geometry
  RIV_cnt <- sum(match(unlist(strsplit(name, ' ')),
                       c(features_RIV_en, features_RIV_fr)), na.rm=T)

  LAC_cnt <- sum(match(unlist(strsplit(name, ' ')),
                       c(features_LAC_en, features_LAC_fr)), na.rm=T)

  if (RIV_cnt>0 & LAC_cnt==0){
    geometry <- "LINE"
  }else if (RIV_cnt==0 & LAC_cnt>0){
    geometry <- "POLYGON"
  }else{
    geometry <- "UNKNOWN"
  }
  return(geometry)
}

