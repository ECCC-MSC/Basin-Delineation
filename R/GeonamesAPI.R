#' Construct search query for CGNDB API
#'
#' @param terms a list() of named search terms, see details for more info
#' @param format Character string; one of c('json', 'html', 'csv', 'kml', 'gml', 'wkt').
#' @details More details about search terms can be found at
#' http://www.nrcan.gc.ca/earth-sciences/geography/place-names/tools-applications/9249
#' Important ones include q, bbox, lat, lon, radius, category, theme, province
#' Theme codes:
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
#' @return Character string; a URL-encoded search string
#' @examples
#' QueryGeogratisAPI()
#' QueryGeogratisAPI(list(q='Ottawa', theme=985), format='csv')
CGNDBQueryAPI <- function(terms, format=NULL){
  args <- lapply(names(terms), function(x) paste(x,'=',terms[x], sep=''))
  args <- paste(unlist(args), collapse='&')
  base.url <- "http://geogratis.gc.ca/services/geoname/en/geonames"
  if (!is.null(format)){
    base.url <- paste(base.url, '.', sep='')
  }
  query <- paste(base.url, format, "?", args, sep="")
  return(query)
}

#' Download data from CGNDB
#'
#' @inheritParams CGNDBQueryAPI
#' @return file path to downloaded file
#' Canadian Geographic Names Database
CGNDBDownloadAPI <- function(dstfile, terms, format){
  url <- CGNDBQueryAPI(terms=terms, format=format)
  download.file(url = url, destfile = dstfile, quiet = FALSE, mode =  "wb", cacheOK = TRUE)
  return(dstfile)
}


# CONCISE CODES: "RIV", "LAKE", "SEA", "CHAN", "BAY",
CGNDBHydroKML <- function(station, dstfile){
  if (!class(station)=='character'){
    stn_name <- station@data$station_name
  }else{
    stn_name <- station
  }
  name <- ParseStationName(stn_name)
  geom <- StationNameGeometry(name)
  if (geom=="LINE"){
    concise="RIV"
  }else if(geom=="POLYGON"){
    concise="LAKE"
  }
  file <- CGNDBDownloadAPI(dstfile=dstfile, terms=list(q=name, theme=979, concise=concise), format='kml')
  file <- FilterKMLGeometry(file, tolower(geom))

  if (is.null(file)){return(NULL)}
  shape <- readOGR(file)


  return(shape)

  # tryCatch({
  # shape <- readOGR(file)
  # return(shape)
  # }, warning = function(w){
  # print(w)
  # return(shape)
  # }, error = function(e){
  # print(e)
  # })
}

FilterKMLGeometry <- function(kml.file, geom, tryAlternate=T){
  geom <- tolower(geom)

  # make sure a feature of the desired type exists...
  feat.kml <- suppressWarnings(readLines(kml.file))
  geom.kml <- switch(geom,
                     'point' ='<Point>',
                     'line'="<LineString>",
                     'polygon'="<Polygon>")
  n.geom <- sum(grepl(geom.kml, feat.kml))  # how many of desired type

  if ((n.geom)==0){
    if (tryAlternate){
      n.poly <- sum(grepl("<Polygon>", feat.kml))
      n.line <- sum(grepl("<Line>", feat.kml))
      n.point <- sum(grepl("<Point>", feat.kml))
    G <- c(n.poly, n.line, n.point)
    if (all(G==0)){return(NULL)}
    geom <- c('polygon', 'line', 'point')[(which(G != 0))[1]]
    }else{
    return(NULL)
  }
  }
  OGR.geom <- switch(geom,
                     'point' ='MultiPoint',
                     'line'="MultiLineString",
                     'polygon'="MultiPolygon")
  out.file <- gsub("\\.kml$", paste("_", tolower(geom), "\\.kml", sep=''), kml.file)

  # convert to desired type
  ogr2ogr(src=kml.file,
          dst=out.file, f='kml',
          where=sprintf("OGR_GEOMETRY = '%s'", OGR.geom),
          mapFieldType = 'DateTime=String')

  return(out.file)
}
