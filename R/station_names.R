
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
  station_name <- toupper(station_name)

  join_terms_en <- c("AT", "ABOVE", "BELOW", "NEAR",
                     "EAST OF", "IN", "OUTFLOW", "BETWEEN", "WEST OF",
                     "NORTH OF", "SOUTH OF")
  join_terms_fr_1 <- c("EN AVAL", "EN AMONT", "CENTRALE", "PRES" )
  join_terms_fr_2 <- c("DU", "DE LA", "DE", "DE L'")
  join_terms_fr_3 <- c("AU", "A LA")

  # try to split on english join term
  splitters <- sapply(join_terms_en, pad, pads = " ")

  # find all splitting terms
  matches <- sapply(splitters,
                    function(p) gregexpr(p, text = station_name)[[1]][1])

  # if any there are any matches, return text before spliting term
  if (!all(matches == -1)){
    matches[matches == -1] <- NA
    first_split <- splitters[which.min(matches)] # first instance of splitting term
    name <- strsplit(station_name, split = first_split)[[1]][1]

  # if no english splitting terms, assume its a french name and move parenthetical description
  }else{
    name <- gsub(
      "([[:alnum:][:blank:]\\.\\'-]*) \\(([[:alnum:][:blank:]\\']*)\\).*",
      "\\2 \\1",
      station_name)
  }

  # Clean up name
  # remove anything after{ a comma, SITE, NO.
  name <- gsub("\\.", "", name)  # strip off any periods
  name <- gsub("' ", "'", name)  # get rid of spaces after apostrophes

  # deal with any remaining parenthetical material?

  return(name)
}

#===============================================================================
#' @title Geometry of waterbody from station name
#'
#' @description Attempts to figure out the geometry of a waterbody (point, polygon,
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

  features_RIV_en <- c("RIVER", "CREEK", "BROOK", "STREAM", "CHANNEL", "CANAL",
                       "RUN", "SPILLWAY", "TRIBUTARY", "DIVERSION", "FLOW")

  features_LAC_fr <- c( "LAC", "BASSIN", "ETANG", "BAIE", "LACS")
  features_LAC_en <- c( "LAKE", "BASIN", "RESERVOIR", "POND", "SWAMP", "SLOUGH",
                        "SWAMP", "BAY", "POND", "LAKES" )

  features_OTH_en <- c("DRAIN", "DITCH", "RILL", "SPRING", "SPRINGS", "COULEE")

  name <- toupper(name)
  name <- gsub("[^[:alnum:] ]", " ", name) # remove wildcards if present

  # Try to determine expected shape geometry
  RIV_cnt <- sum(match(unlist(strsplit(name, " ")),
                       c(features_RIV_en, features_RIV_fr)), na.rm = T)

  LAC_cnt <- sum(match(unlist(strsplit(name, " ")),
                       c(features_LAC_en, features_LAC_fr)), na.rm = T)

  if (RIV_cnt > 0 & LAC_cnt == 0){
    geometry <- "LINE"
  }else if (RIV_cnt == 0 & LAC_cnt > 0){
    geometry <- "POLYGON"
  }else{
    geometry <- "UNKNOWN"
  }

  return(geometry)
}


#==============================================================================
#' @title Add character strings on either end of an input character string
#'
#' @param text input character string
#'
#' @param pads character string to append and prepend to text
#'
#' @param rpads (optional) if provided, adds different text on right-hand side
#' of string.
#'
#' @return character string, original text with pads appended and prepended
#'
#' @examples
#' pad('file', '00')
#' pad('filename', '', '.csv')
#'
#' @export
#'
#' @keywords internal
#==============================================================================
pad <- function(text, pads, rpads){
  if (missing(rpads)){
    sprintf(paste("%s", text, "%s", sep = ""), pads, pads)
  }else{
    sprintf(paste("%s", text, "%s", sep = ""), pads, rpads)
  }
}
