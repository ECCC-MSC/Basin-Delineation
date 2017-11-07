#con <- dbConnect(RSQLite::SQLite(), "M:\\trasnfer\\Hydat_august\\Hydat.sqlite3")

#' Length of Longitudinal Degree
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

#' Which UTM Zone am I in?
#'
#' @description Determines UTM zone from longitude
#' @param long numeric longitude in decimal degrees
#' @return integer UTM zone number
#' @export
WhichZone <- function(long){
  zone <- (floor((long + 180)/6) %% 60) + 1
  return(zone)
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

#' HydroSHEDS DEM index
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

#' HydroMosaic
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


#' Coordinate systems
#'
#' @description A helper function to easily produce proj4 strings without introducting global variables
#' and without cluttering the main body of the code with long proj4 strings.
#' @param x character string of a
#' @details
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


TileIndex <- function(geom1, tile, tilename){
  # read-in tile if filename is provided
  tile <- InterpretShapefile(tile)
  # check CRS, transform if necessary
  geom1 <- SameCRS(geom1,tile)
  # get intersecting tile iDs
  tiles <- over(geom1, tile, returnList = T)[[1]]
  tiles <- tiles[,tilename]
  return(as.character(tiles))
}

#' Check for same CRS between Spatial* objects
#'
#' @description Checks whether or not two objects have the same spatial reference.  If not, one of them
#' is transformed to match the other.
#' @param spgeom1 Spatial* object that will be evaluated and transformed if necessary
#' @param spgeom2 Spatial* object that will not be transformed
SameCRS <- function(spgeom1, spgeom2){
  if (spgeom1@proj4string@projargs != spgeom2@proj4string@projargs){
    spgeom1 <- sp::spTransform(spgeom1, CRSobj = sp::CRS(spgeom2@proj4string@projargs))
  }
  return(spgeom1)
}

#' Shapefile Helper
#'
#' @description allows for shapefiles to be passed as character strings or as R spatial objects in other functions
#' @param x either an R spatial object or a character string specifying a shapefile path
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

#' Bounding box buffer
#'
#' @description Increases the size of a bounding box by a specified
#' @param geom1 an R Spatial* object
#' @param tol numeric, distance in kilometers to buffer
#' @return a bounding box R object
#' @details This is used with \code{link[rcanvec]{nts.bbox}} in order to get all NTS tiles that
#' are within the buffer distance of the object
ExpandBBox <- function(geom1, tol){
  geom1 <- sp::spTransform(geom1, GetProj4("WGS84"))
  box <- sp::bbox(geom1)
  dlat <- (tol * 1e-3) / 111
  dlon <- (tol * 1e-3) / LongitudeLength(mean(box["latitude",]))
  box.new <- box + c(-dlon, -dlat, dlon, dlat)
  return(box.new)
}


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
