#con <- dbConnect(RSQLite::SQLite(), "M:\\trasnfer\\Hydat_august\\Hydat.sqlite3")

#' Calculate Length of Longitudinal Degree
#'
#' @description  Calculates the length of one degree of longitude at a given latitude
#' @param lat  numeric vector - northing of location (degrees or decimal degrees)
#' @param lat.m (optional) numeric vector - northing of location, minutes (if lat is not specified in decimal degrees)
#' @param lat.s (optional) numeric vector - northing of location, seconds (if lat is not specified in decimal degrees)
#' @return numeric vector of lengths  (measured in kilometres)
#' @export
LongitudeLength <- function(lat, lat.m, lat.s){
  lon.at.eq = 111.321
  mm <- 0
  ss <- 0
  if (!missing(lat.m)){
    mm <- lat.m/60
  }
  if (!missing(lat.s)){
    ss <- lat.s/3600
  }
  lat.dd <- lat + mm + ss
  lat.rad <- lat.dd * pi/180

  lon.km = cos(lat.rad) * lon.at.eq
  return(lon.km)
}

#' Calculate UTM zone from longitude
#'
#' @description Determines UTM zone from longitude
#' @param long numeric longitude in decimal degrees
#' @details (floor((long + 180)/6) %% 60) + 1
#' @return integer UTM zone number
#' @export
WhichZone <- function(long){
  zone <- (floor((long + 180)/6) %% 60) + 1
  return(zone)
}


#' Get index of HydroSHEDS DEM tile
#'
#' @description returns the name of a HYDROSHEDS tiled data file name based on the coordinates of a point
#' @param long numeric longitude
#' @param lat numeric latitude
#' @param fext file extension for output filename (needs leading '.')
#' @keywords internal
#' @export
HydroTile <- function(long, lat, fext=''){
  x <- floor(long/5)*5
  y <- floor(lat/5)*5
  xd <- ifelse(x<0, 'w', 'e')
  yd <- ifelse(y<0, 's', 'n')
  out <- sprintf("%s%02.2d%s%03.3d_con%s",yd, abs(y),xd, abs(x), fext)
  return(out)
}

#' Find all hydroSHEDS tiles within a distance of a point
#'
#' @description Finds all HydroSHEDS DEM tiles within a specified tolerance of a point.  To be used
#' with 3 arc-second hydrologically conditioned DEM.
#' @param tol Tolerance in m to consider
#' @return A vector of grid names to load
#' @keywords internal
#' @export
HydroMosaic <- function(long, lat, tol, ...){
  tol <- tol * 0.001 # convert to km
  LL <- LongitudeLength(lat)
  grids <- HydroTile(long, lat, ...)

  # check all sides for proximity, add grids if necessary
  if ((lat - floor(lat / 5) * 5) * 111 < tol){ # too close to bottom
    grids <- c(grids, HydroTile(long, lat - 5, ...))
  }else if ((ceiling(lat / 5) * 5 - lat) * 111 < tol){ # too close to top
    grids <- c(grids, HydroTile(long, lat + 5, ...))
  }
  if ((long-floor(long / 5) * 5) * LL < tol){ # too far west
    grids <- c(grids, HydroTile(long - 5, lat, ...))
  }else if ((ceiling(long / 5) * 5 - long)*LL < tol){
    grids <- c(grids, HydroTile(long + 5, lat, ...))
  }
  return(grids)
}

#' Get Mosaic Extent
#'
#' @param names names of mosaic file
GetMosaicLimits <- function(names){
  names <- gsub("_con.*", "", names)
  names <- gsub("[ns]", "", names)
  names <- gsub("[ew]", ",", names)
  names <- strsplit(names,",")
  return(names)
}

#' Read SAGA Grid header
#'
#' @description get info for SAGA grid
#' @param file character string path to *.sgrd
#' @return data frame with list of grid parameters
#' @export
#' @keywords internal
read.sgrd.header <- function(file){
  # get .sgrd grid info from ascii header
  data <- read.table(file, sep=c('='), strip.white = T, stringsAsFactors = F)
  output <- data.frame(t(data[,2]), stringsAsFactors = F)
  names(output) = data[,1]
  #output <- lapply(output, type.convert)
  return(output)
}


#' Proj4 string from abbreviation
#'
#' @description A helper function to easily produce proj4 strings without introducting global variables
#' and without cluttering the main body of the code with long proj4 strings.
#' @param x character string of a coordinate system. See options in details below
#' @return proj4 string corresponding to x
#' @keywords internal
#' @details
#'  "WGS84" = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
#'  "AlbersEqualAreaConic" = "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40
#'  +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#'  ... more to come
GetProj4 <- function(x){
  return(
    switch(x,
      "WGS84"="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
      "AlbersEqualAreaConic"="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
    )
  )
}

#'#' Hybas to stationprefix
#'
#' @description For a station prefix, returns the first 3 digits of the pfaffstetter code corresponding
#' to the nearest hydrobasins polygon. This allows the hydrobasins file to be subset using indexing before
#' running through searches by removing polygons that are neither upstream nor downstream of the station
#' @param station character string specifying a station name e.g. "02AB001"
#' @keywords internal
NearestHYBAS <- function(station){
  code = regmatches(x =station ,m = regexpr("^\\d{2}", station))
  basins <- switch(code,
             "01"=c(723, 726),
             "02"=c(725, 724, 723, 726),
             "03"=c(721, 722, 723, 724),
             "04"=c(713),
             "05"=c(712),
             "06"=c(711, 832, 833, 852, 851),
             "07"=c(822),
             "08"=c(783, 782),
             "09"=c(811, 812),
             "10"=c(812, 822, 823, 831, 841, 842, 843, 844, 861, 862, 863),
             "11"=c(712, 742))
  return(basins)
}

GetTilePathsHS <- function(names, DEM.dir){
  original.DEM <- unlist(lapply(names, list.files, path=DEM.dir, full.names = T, recursive = T, include.dirs = T))
  original.DEM <- original.DEM[basename(original.DEM)==basename(dirname(original.DEM))]
  return(original.DEM)
}

#' Get names of DEM tiles that overlap object
#' @description Find all DEM tiles that intersect a spatial object
#' @param geom1 An R Spatial* object
#' @param tileindex either an R SpatialPolygonsDataFrame of the DEM tile index, or a character
#' path pointing to such a shapefile
#' @param tile.id.field Name of column in tileindex that gives DEM sheet number
#' @export
#' @keywords internal
TileIndex <- function(geom1, tileindex, tile.id.field){
  # read-in tile if filename is provided
  tile <- InterpretShapefile(tileindex)
  # check CRS, transform if necessary
  geom1 <- SameCRS(geom1,tile)
  # get intersecting tile iDs
  tiles <- over(geom1, tile, returnList = T)[[1]]
  tiles <- tiles[,tile.id.field]
  return(as.character(tiles))
}

#' Ensure that two Spatial* objects have the same CRS
#'
#' @description Checks whether or not two objects have the same spatial reference.  If not, one of them
#' is transformed to match the other.
#' @param spgeom1 Spatial* object that will be evaluated and transformed if necessary
#' @param spgeom2 Spatial* object that will not be transformed
#' @export
SameCRS <- function(spgeom1, spgeom2){
  if (spgeom1@proj4string@projargs != spgeom2@proj4string@projargs){
    spgeom1 <- sp::spTransform(spgeom1, CRSobj = sp::CRS(spgeom2@proj4string@projargs))
  }
  return(spgeom1)
}

#' Shapefile Helper
#'
#' @description allows for shapefiles to be passed as character strings or
#' as R spatial objects in other functions
#' @param x either an R spatial object or a character string specifying a shapefile path
#' @export
#' @keywords internal
InterpretShapefile <- function(x){
  if (class(x)=="character"){
    x <- rgdal::readOGR(x)
    return(x)
  }
  if (grepl("Spatial", class(x))){
    return(x)
  }else{
    stop("Could not interpret shapefile")
  }
}

#' Expand Bounding box
#'
#' @description Increases the size of a bounding box by a specified
#' @param geom1 an R Spatial* object
#' @param tol numeric, distance in kilometers to buffer
#' @return a bounding box R object
#' @details This is used with \code{link[rcanvec]{nts.bbox}} in order to get all NTS tiles that
#' are within the buffer distance of the object
#' @export
#' @keywords internal
ExpandBBox <- function(geom1, tol){
  geom1 <- sp::spTransform(geom1, GetProj4("WGS84"))
  box <- sp::bbox(geom1)
  dlat <- (tol * 1e-3) / 111
  dlon <- (tol * 1e-3) / LongitudeLength(mean(box[2,]))
  if (!(all(0 < abs(box[2,])) & all(abs(box[2,]) < 90))){stop()}
  box.new <- box + c(-dlon, -dlat, dlon, dlat)
  return(box.new)
}

#' Snap point to nearest line
#' @description moves a point to the nearest point on a line
#' @param point spatialPointsDataFrame
#' @param lines spatialLinesDataFrame or spatialLines
#' @return a spatial point with the same data as the input point, but on one of the lines in lines.
#' @export
#' @keywords internal
SnapToNearest <- function(point, lines){
  point <- SameCRS(point, lines) # put them in same CRS
  n <- rgeos::gNearestPoints(point, lines)[2]  #first element is original point
  pointNew <- SpatialPointsDataFrame(coords=n@coords, data=point@data, proj4string = point@proj4string)
  return(pointNew)
}

#' Create SpatialPointsDataFrame from csv
#'
#' @description Creates an R SpatialPointsDataFrame from a csv or
#' @param file path to *.csv file
#' @param headerX character, column name containing x coordinate information
#' @param headerY character, column name containing y coordinate information
#' @param CRS (optional) character, proj4 string specifying the projection information of the points. If missing,
#' the default is WGS84
#' @return SpatialPointsDataFrame
#' @export
SpatialCSV <- function(file, headerX, headerY, CRS){
  if (missing(CRS)){
    CRS <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  }
  CRS <- sp::CRS(CRS)
  file.data <- read.csv(file)
  coords <- file.data[,c(headerX, headerY)]
  output <- SpatialPointsDataFrame(coords, data=file.data, proj4string = CRS)
  return(output)
}


#' Find degree-graticule tiles within distance
#'
GraticuleIndices <- function(long, lat, tol, resolution=1, ...){
  tol <- tol * 0.001 # convert to km
  LL <- LongitudeLength(lat)
  dlat <- tol/111  #convert to degree
  dlon <- tol/LL   #convert to degree

  xrange <- c(long-dlon, long+dlon)
  xrange <- floor(xrange/resolution)*resolution
  xrange <- seq(min(xrange), max(xrange), resolution)
  yrange <- c(lat-dlat, lat+dlat)
  yrange <- ceiling(yrange/resolution)*resolution
  yrange <- seq(min(yrange), max(yrange), resolution)

  out <- expand.grid(xrange, yrange)
  return(out)
}

#only works with points, and only works with lat/lon
NEDcoverage <- function(geom1, tol, ...){
  if (grepl("units=m", geom1@proj4string@projargs)){
    geom1 <- sp::spTransform(geom1,
                                CRSobj = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  }
  coords <- gCentroid(geom1)@coords[1,]
  ind <- GraticuleIndices(coords[1], coords[2], tol=tol, ...)
  apply(ind, 1, function(x) USGSTileName(x[1], x[2]))
}
