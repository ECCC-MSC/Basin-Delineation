#===============================================================================
#' @title Construct query for CGNDB API
#'
#' @description Returns a URL that can be used to query the
#' Canadian Geographic Names Database using the online API
#'
#' @param terms a list of named search terms, see details for more info
#'
#' @param format Character string; one of c('json', 'html', 'csv', 'kml', 'gml',
#'  'wkt'). Determines the format of the results returned by the API.
#'
#' @details The terms parameter must be a list with the names corresponding to
#' search parameters of the CGNDB.
#' More details about search terms can be found at
#' http://www.nrcan.gc.ca/earth-sciences/geography/place-names/tools-applications/9249
#' Common codes include q, bbox, lat, lon, radius, category, theme, and province
#'
#' Theme codes allow for only a certain category of feature to be returned. The
#' theme codes are:
#' 979  Water Feature
#' 980 	Volcanic Feature
#' 981 	Undersea and Maritime Feature
#' 982 	Underground Feature
#' 983 	Unclassified
#' 984 	Terrain Feature
#' 985 	Populated Place
#' 986 	Ice and Snow Feature
#' 987 	Feature Associated with Vegetation
#' 988 	Constructed Feature
#' 989 	Administrative Area
#'
#' @return Character string; a URL-encoded search string
#'
#' @examples
#' QueryGeogratisAPI()
#' QueryGeogratisAPI(list(q='Ottawa', theme=985), format='csv')
#===============================================================================
CGNDBQueryAPI <- function(terms, format=NULL){
  args <- lapply(names(terms), function(x) paste(x, "=", terms[x], sep = ""))
  args <- paste(unlist(args), collapse = "&")
  base.url <- "http://geogratis.gc.ca/services/geoname/en/geonames"
  if (!is.null(format)){
    base.url <- paste(base.url, ".", sep = "")
  }
  query <- paste(base.url, format, "?", args, sep = "")
  return(query)
}

#==============================================================================
#' @title Download data from CGNDB
#'
#' @inheritParams CGNDBQueryAPI
#'
#' @param dstfile character string specifying the path to a destination file
#' to store query results
#'
#' @description Generates a query and downloads the results
#'
#' @return file path to downloaded file from Canadian Geographic Names Database
#'
#' @export
#==============================================================================
CGNDBDownloadAPI <- function(dstfile, terms, format){
  url <- CGNDBQueryAPI(terms = terms, format = format)
  download.file(url = url, destfile = dstfile, quiet = FALSE, mode =  "wb",
                cacheOK = TRUE)
  return(dstfile)
}

#===============================================================================
#' @title download kml of hydrological feature
#'
#' @description Based on the name of a NHS station, attempts to download a
#' kml corresponding to the hydrological feature being measured.
#'
#' @details After downloading, the kml is converted to a shapefile. Because of the
#' limitations of feature types in a shapefile, only one geometry type is preserved
#' (points, lines or polygons). The script does its best to figure out
#' which of the three should be preserved based on the station name.
#'
#' @param station name of a NHS / HYDAT hydrometric station or the name of a
#' waterbody
#'
#' @param dstfile character string specifying a kml file to download
#'
#' @details CONCISE CODES: "RIV", "LAKE", "SEA", "CHAN", "BAY",
#'
#' @return a Spatial* object
#===============================================================================
CGNDBHydroKML <- function(station, dstfile){

  if (!class(station) == "character"){
    stn_name <- station@data$station_name
  }else{
    stn_name <- station
  }

  name <- ParseStationName(stn_name)
  geom <- StationNameGeometry(name)

  if (geom == "LINE"){
    concise <- "RIV"
  }else if (geom == "POLYGON"){
    concise <- "LAKE"
  }
  file <- CGNDBDownloadAPI(dstfile = dstfile,
                           terms = list(q = name, theme = 979,
                                          concise = concise), format = "kml")
  file <- FilterKMLGeometry(file, tolower(geom))

  if (is.null(file)){ return(NULL) }

  shape <- rgdal::readOGR(file)

  return(shape)
}

#===============================================================================
#' @title Convert KML to shapefile, retaining only one geometry type
#'
#' @param kml.file character path to kml file
#'
#' @param geom character string, which feature type do you want to retain? options
#' are one of ('point', 'line', 'polygon').
#'
#' @param try_alternate logical, if no features of the desired type are found,
#' whether or not to try loading a different feature type.
#'
#' @return character, path to newly created shapefile
#===============================================================================
FilterKMLGeometry <- function(kml.file, geom, try_alternate=T){
  geom <- tolower(geom)

  # make sure a feature of the desired type exists...
  feature_kml <- suppressWarnings(readLines(kml.file))

  kml_geom <- switch(geom,
                     "point" = "<Point>",
                     "line" = "<LineString>",
                     "polygon" = "<Polygon>")

  n_geom <- sum(grepl(kml_geom, feature_kml))  # how many of desired type

  if ( (n_geom) == 0){
    if (try_alternate){
      n_poly <- sum(grepl("<Polygon>", feature_kml))
      n_line <- sum(grepl("<Line>", feature_kml))
      n_point <- sum(grepl("<Point>", feature_kml))
    G <- c(n_poly, n_line, n_point)

    if (all(G == 0)){ return(NULL) }

    geom <- c("polygon", "line", "point")[(which(G != 0))[1]]
    }else{
    return(NULL)
      }
  }

  OGR_geom <- switch(geom,
                     "point" = "MultiPoint",
                     "line" = "MultiLineString",
                     "polygon" = "MultiPolygon")
  out_file <- gsub("\\.kml$", paste("_", tolower(geom), "\\.shp", sep = ""),
                   kml.file)

  # convert to desired type
  gdalUtils::ogr2ogr(src = kml.file,
          dst = out_file, f = "ESRI Shapefile",
          where = sprintf("OGR_GEOMETRY = '%s'", OGR_geom),
          mapFieldType = "DateTime=String")

  return(out_file)
}
