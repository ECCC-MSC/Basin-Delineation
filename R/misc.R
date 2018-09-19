#==============================================================================
#' @title Calculate Length of Longitudinal Degree
#'
#' @description  Calculates the length of one degree of longitude at a given latitude
#'
#' @param lat  numeric vector - northing of location (degrees or decimal degrees)
#'
#' @param lat_m (optional) numeric vector - northing of location, minutes (if lat is not specified in decimal degrees)
#'
#' @param lat_s (optional) numeric vector - northing of location, seconds (if lat is not specified in decimal degrees)
#'
#' @return numeric vector of lengths  (measured in kilometres)
#'
#' @examples
#' LongitudeLength(61.4)
#' LongitudeLength(34,2,24)
#'
#' @export
#==============================================================================
LongitudeLength <- function(lat, lat_m, lat_s){
  lon_at_eq <- 111.321
  mm <- 0
  ss <- 0
  if (!missing(lat_m)){
    mm <- lat_m / 60
  }
  if (!missing(lat_s)){
    ss <- lat_s / 3600
  }
  lat_dd <- lat + mm + ss
  lat_rad <- lat_dd * pi / 180

  lon_km <- cos(lat_rad) * lon_at_eq
  return(lon_km)
}

#==============================================================================
#' @title Determine UTM zone from longitude
#'
#' @description Determines UTM zone from longitude
#'
#' @param long numeric longitude in decimal degrees
#'
#' @details (floor((long + 180)/6) %% 60) + 1
#'
#' @return integer UTM zone number
#'
#' @examples
#' WhichZone(-125.7)
#'
#' @export
#==============================================================================
WhichZone <- function(long){
  zone <- (floor( (long + 180) / 6) %% 60) + 1
  return(zone)
}

#==============================================================================
#' @title Get name of HydroSHEDS DEM tile
#'
#'
#' @description returns the name of a HYDROSHEDS tiled data file name based on
#' the coordinates of a point
#'
#' @param long numeric longitude
#'
#' @param lat numeric latitude
#'
#' @param fext file extension for output filename
#'
#' @keywords internal
#'
#' @examples
#' HydroTile(-126.4, 56)
#'
#' @export
#==============================================================================
HydroTile <- function(long, lat, fext=""){
  fext <- gsub("^([^\\.])", "\\.\\1", fext)

  x <- floor(long / 5) * 5
  y <- floor(lat / 5) * 5
  if (y >= 60){stop("No HydroSHEDS data north of 60 degrees")}
  xd <- ifelse(x < 0, "w", "e")
  yd <- ifelse(y < 0, "s", "n")
  out <- sprintf("%s%02.2d%s%03.3d_con%s", yd, abs(y), xd, abs(x), fext)

  return(out)
}

#==============================================================================
#' Find all hydroSHEDS tiles within a distance of a point
#'
#' @description Finds all HydroSHEDS DEM tiles within a specified tolerance of
#' a point.  To be used with 3 arc-second hydrologically conditioned DEM.
#' @param long longitude in decimal degrees
#' @param lat latitude in decimal degrees
#' @param tol Tolerance in m to consider
#' @return A vector of grid names to load
#' @keywords internal
#' @examples
#' HydroMosaic(-115, 52, 50000)
#' @export
#==============================================================================
# HydroMosaic <- function(long, lat, tol, ...){
#   tol <- tol * 0.001 # convert tolerance to km
#   LL <- LongitudeLength(lat)
#   grids <- HydroTile(long, lat, ...)
#
#   # check all sides for proximity, add grids if necessary
#
#   if ( (lat - floor(lat / 5) * 5) * 111 < tol){
#     # too close to bottom
#     grids <- c(grids, HydroTile(long, lat - 5, ...))
#   }else if ( (ceiling(lat / 5) * 5 - lat) * 111 < tol){
#     # too close to top
#     grids <- c(grids, HydroTile(long, lat + 5, ...))
#   }
#
#   if ( (long - floor(long / 5) * 5) * LL < tol){
#     # too far west
#     grids <- c(grids, HydroTile(long - 5, lat, ...))
#   }else if ( (ceiling(long / 5) * 5 - long) * LL < tol){
#     grids <- c(grids, HydroTile(long + 5, lat, ...))
#   }
#
#   return(grids)
# }

HydroMosaic <- function(long, lat, tol, ...){
  tol   <- tol * 0.001 # convert tolerance to km
  # convert to degrees
  toly  <- tol / 111
  tolx  <- tol / LongitudeLength(lat)

  # check all sides for proximity, add grids if necessary
  xmin <- floor((long - tolx) / 5) * 5
  xmax <- floor((long + tolx) / 5) * 5

  xrange <- seq(xmin, xmax, 5)
  yrange <- seq(floor((lat - toly) / 5) * 5, floor((lat + toly) / 5) * 5, 5)

  pairs <- expand.grid(xrange, yrange)

  grids <- apply(pairs, 1, function(x) HydroTile(x[1], x[2]))
  print(grids)
  return(grids)
}


#' Get Mosaic Extent
#'
#' @param names names of mosaic file
GetMosaicLimits <- function(names){

  names <- gsub("_con.*", "",  names)
  names <- gsub("[ns]",   "",  names)
  names <- gsub("[ew]",   ",", names)

  names <- strsplit(names, ",")

  return(names)
}

#==============================================================================
#' @title Read SAGA SGRD header
#'
#' @description get info for SAGA grid
#'
#' @param file character string path to *.sgrd
#'
#' @return data frame with list of grid parameters
#'
#' @export
#'
#' @keywords internal
#==============================================================================
read.sgrd.header <- function(file){
  # get .sgrd grid info from ascii header
  data <- read.table(file, sep = c("="), strip.white = T, stringsAsFactors = F)
  output <- data.frame(t(data[, 2]), stringsAsFactors = F)
  names(output) <- data[, 1]

  return(output)
}

#==============================================================================
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
#==============================================================================
GetProj4 <- function(x){
  return(
    switch(x,
      "WGS84" = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
      "AlbersEqualAreaConic" = "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
    )
  )
}

#==============================================================================
#'Hybas to stationprefix
#'
#' @description For a station prefix, returns the first 3 digits of the pfaffstetter code corresponding
#' to the nearest hydrobasins polygon. This allows the hydrobasins file to be subset using indexing before
#' running through searches by removing polygons that are neither upstream nor downstream of the station
#'
#' @param station character string specifying a station name e.g. "02AB001"
#'
#' @keywords internal
#==============================================================================
NearestHYBAS <- function(station){
  if (!grepl("^\\d{2}[[:alpha:]]{2}\\d{3}", station)){
    return(seq(100, 999))
  }
  code <- regmatches(x = station, m = regexpr("^\\d{2}", station))
  basins <- switch(code,
             "01" = c(723, 726),
             "02" = c(725, 724, 723, 726),
             "03" = c(721, 722, 723, 724),
             "04" = c(713),
             "05" = c(712),
             "06" = c(711, 832, 833, 852, 851),
             "07" = c(822),
             "08" = c(783, 782),
             "09" = c(811, 812),
             "10" = c(812, 822, 823, 831, 841, 842, 843, 844, 861, 862, 863),
             "11" = c(712, 742))
  return(basins)
}

#==============================================================================
#' @title Get file paths of hydrosheds tiles
#'
#'
#' @param names names of hydrosheds tiles (e.g. 'n45w060')
#'
#' @param DEM_dir character path pointing to directory containing hydrosheds
#' tiles
#'
#' @return character, list of file paths for each of the tiles specified
#' in the names parameter
GetTilePathsHS <- function(names, DEM_dir){
  original_DEM <- unlist(lapply(
    names, list.files, path = DEM_dir, full.names = T, recursive = T,
    include.dirs = T))

  original_DEM <- original_DEM[basename(original_DEM) ==
                                 basename(dirname(original_DEM))]
  return(original_DEM)
}

#==============================================================================
#' @title Get names of DEM tiles that overlap object
#'
#' @description Find all DEM tiles that intersect a spatial object
#'
#' @param geom1 An R Spatial* object
#'
#' @param tileindex either an R SpatialPolygonsDataFrame of the DEM tile index, or a character
#' path pointing to such a shapefile
#'
#' @param tileid_field Name of column in tileindex that gives DEM sheet number
#'
#' @export
#'
#' @keywords internal
#==============================================================================
TileIndex <- function(geom1, tileindex, tileid_field){
  # read-in tile if filename is provided
  tile <- InterpretShapefile(tileindex)

  # check CRS, transform if necessary
  geom1 <- SameCRS(geom1, tile)

  # get intersecting tile iDs
  tiles <- sp::over(geom1, tile, returnList = T)[[1]]
  tiles <- tiles[, tileid_field]
  return(as.character(tiles))
}

#==============================================================================
#' @title Ensure that two Spatial* objects have the same CRS
#'
#' @description Checks whether or not two objects have the same spatial
#'  reference.  If not, one of them is transformed to match the other.
#'
#' @param spgeom1 Spatial* object that will be evaluated and transformed
#'  if necessary
#'
#' @param spgeom2 Spatial* object that will not be transformed (reference CRS)
#'
#' @return a Spatial* object with the same CRS as spgeom2
#'
#' @export
#==============================================================================
SameCRS <- function(spgeom1, spgeom2){
  if (spgeom1@proj4string@projargs != spgeom2@proj4string@projargs){
    spgeom1 <- sp::spTransform(spgeom1, CRSobj = sp::CRS(spgeom2@proj4string@projargs))
  }

  return(spgeom1)
}

#==============================================================================
#' @title read shapefile as filepath or R object
#'
#' @description allows for shapefiles to be passed as character strings or
#' as R spatial objects in other functions.
#'
#' @param x either an R spatial object or a character string specifying a
#' shapefile path
#'
#' @param use_sf logical, whether or not to try to use use the 'sf' package
#' to read the shapefile quickly.
#'
#' @export
#'
#' @keywords internal
#==============================================================================
InterpretShapefile <- function(x, use_sf=T, quiet=T){
  if (class(x) == "character"){

    if (require(sf)){
      x <- sf::st_read(x, quiet=quiet)
      x <- sf::st_zm(x)  # drop Z and M dimensions
      x <- as(x, "Spatial")
    }else{
      x <- rgdal::readOGR(x, verbose=!quiet)
    }

    return(x)
  }

  if (grepl("Spatial", class(x))){
    return(x)
  }else{
    stop("Could not interpret shapefile")
  }
}

#==============================================================================
#' Expand Bounding box
#'
#' @description Increases the size of a bounding box by a specified
#'
#' @param geom1 an R Spatial* object
#'
#' @param tol numeric, distance in kilometers to buffer
#'
#' @return a bounding box R object
#'
#' @details This is used with \code{link[rcanvec]{nts.bbox}} in order to get all NTS tiles that
#' are within the buffer distance of the object
#'
#' @export
#'
#' @keywords internal
#==============================================================================
ExpandBBox <- function(geom1, tol){
  geom1 <- sp::spTransform(geom1, GetProj4("WGS84"))
  box <- sp::bbox(geom1)
  dlat <- (tol * 1e-3) / 111
  dlon <- (tol * 1e-3) / LongitudeLength(mean(box[2, ]))
  if ( !(all(0 < abs(box[2, ])) & all(abs(box[2, ]) < 90))){stop()}
  box_new <- box + c(-dlon, -dlat, dlon, dlat)
  return(box_new)
}

#==============================================================================
#' @title Snap point to nearest line
#'
#' @description moves a point to the nearest point on a line
#'
#' @param point spatialPointsDataFrame
#'
#' @param lines spatialLinesDataFrame or spatialLines
#'
#' @return a spatial point with the same data as the input point, but on one of
#' the lines in lines.
#'
#' @export
#'
#' @keywords internal
#==============================================================================
SnapToNearest <- function(point, lines){
  point <- SameCRS(point, lines) # put them in same CRS
  n <- rgeos::gNearestPoints(point, lines)[2]  #first element is original point
  new_point <- SpatialPointsDataFrame(coords = n@coords, data = point@data,
                                     proj4string = point@proj4string)

  return(new_point)
}




#==============================================================================
#' @title Find intersecting graticule tiles within distance of point
#'
#' @param long target longitude
#'
#' @param lat target latitude
#'
#' @param tol distance (in metres) to search. Search area is a square.
#'
#' @param resolution numeric value indicating the width and height of
#' graticule cells. For example a resolution of 1 corresponds to 1x1 degree tiles
#' and a resolution of 5 corresponds to 5x5 degree tiles
#'
#' @details The graticules are reported based on their upperleftmost corner
#'
#' @return a dataframe whose rows correspond to graticule index coordinates
#'
#==============================================================================
GraticuleIndices <- function(long, lat, tol, resolution=1){
  tol <- tol * 0.001 # convert to km

  LL <- LongitudeLength(lat)
  dlat <- tol / 111  #convert to degree
  dlon <- tol / LL   #convert to degree

  xrange <- c(long - dlon, long + dlon)
  xrange <- floor(xrange / resolution) * resolution
  xrange <- seq(min(xrange), max(xrange), resolution)

  yrange <- c(lat - dlat, lat + dlat)
  yrange <- ceiling(yrange / resolution) * resolution
  yrange <- seq(min(yrange), max(yrange), resolution)

  out <- expand.grid(xrange, yrange)

  return(out)
}

#==============================================================================
#only works with lat/lon

#' @param geom1 geometry with which to cover with NED DEM tiles
#'
#' @param tol distance to buffer geometry
#==============================================================================
NEDcoverage <- function(geom1, tol){
  if (grepl("units=m", geom1@proj4string@projargs)){
    geom1 <- sp::spTransform(geom1,
           CRSobj = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  }

    # for points
    if (all(geom1@bbox[, 1] == geom1@bbox[, 2])){
      coords <- rgeos::gCentroid(geom1)@coords[1, ]
      ind <- GraticuleIndices(coords[1], coords[2], tol = tol, resolution = 1)


      }else{
      # for non-points
      xrange <- c(seq(geom1@bbox[1, 1], geom1@bbox[1, 2], 1), geom1@bbox[1, 2]) # duplicate last to ensure its included
      yrange <- c(seq(geom1@bbox[2, 1], geom1@bbox[2, 2], 1), geom1@bbox[2, 2])
      allpts <- expand.grid(unique(xrange), unique(yrange))
      ind <- apply(allpts, 1, function(x) GraticuleIndices(x[1], x[2], tol=tol, resolution = 1))
      ind <- unique(do.call('rbind', ind))
      }

  # use graticule indices to build tile names
  tiles <- apply(ind, 1, function(x) USGSTileName(x[1], x[2]))

  return(as.character(tiles))
}

#==============================================================================
#' @title Add missing columns
#'
#' @param df A dataframe
#'
#' @param columns a character vector of columns that are desired in output
#' dataframe
#'
#' @param return a dataframe containing all columns specified in columns
#' parameter
#==============================================================================
AddMissingColumns <- function(df, columns){
  missing.cols <- !(columns %in% names(df))
  to_add <- as.data.frame(matrix(nrow = nrow(df), ncol = sum(missing.cols)))
  names(to_add) <- columns[missing.cols]
  df <- cbind(df, to_add)
  return(df)
}

#==============================================================================
#' @title Remove holes from polygon
#'
#' @description Returns the largest ring from a single polygon feature,
#' retaining original ID value.
#'
#'  @details Finds the biggest ring of a SpatialPolygons* object and uses that
#' ring to create a new SpatialPolygons object with the same ID value. This can
#'  be used to clean up holes/overlap within polygons, or to remove small
#'  polygons within a multipart polygon
#'
#' @param poly A SpatialPolygons* object. The feature may have multiple rings.
#' If an object with multiple shapes is passed, the function is iterated over
#' all entries.
#'
#' @return SpatialPolygon
#==============================================================================
outerRing <- function(poly){
  reduce <- function(x){
    ring_num <- which.max( # keep biggest ring
    sapply(
      x@polygons[[1]]@Polygons, slot, name = "area"
    ))
  original.ID <- sapply(slot(x, "polygons"), function(x) slot(x, "ID"))
  ring <- SpatialPolygons(
    list(
      Polygons(
        list(
          x@polygons[[1]]@Polygons[[ring_num]]
        ),
        ID = original.ID)), proj4string = x@proj4string)
  return(ring)
  }

  if (length(poly@polygons) == 1){
    return(reduce(poly))
  }else if (length(poly@polygons) > 1){
    L <- lapply(seq_along(poly), function(S) reduce(poly[S, ]))
    L <- do.call(rbind, L)
    return(L)
    }
  }

zeroBuffer <- function(SpatialObj, byid=TRUE){

  if ("data" %in% slotNames(SpatialObj)){
    fixed_geom <-rgeos::gBuffer(SpatialObj, width = 0, byid = byid)

    f <- switch(class(SpatialObj),
                "SpatialPolygonsDataFrame" = sp::SpatialPolygonsDataFrame,
                "SpatialPointsDataFrame" = sp::SpatialPointsDataFrame,
                "SpatialLinesDataFrame" = sp::SpatialLinesDataFrame)
    out <- f(Sr = fixed_geom, data = SpatialObj@data)

    }else{
      out <- datargeos::gBuffer(SpatialObj, width = 0, byid = byid)
    }
   return(out)
}


#==============================================================================
#' @title Fill gaps between catchment and HydroBASIN
#'
#' @param station spatialpoint*
#'
#' @param basin spatialpolygon* of contributing area basin, already clipped to hybas
#'
#' @param hybas spatialpolygon* corresponding to the entire hydrobasins boundary
#'
#' @param additive whether to add gaps to basin or subtract the downstream gap
#' from the containing hydrobasin
#'
#' @export
#==============================================================================
fill_upstream_gaps <- function(station, basin, hybas, additive=T,
                               cl_strategy="POLYGONATION"){
  print('filling upstream gaps')
  albers <- sp::CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0
                           +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  bprj <- basin@proj4string

  # Ensure consistent projected CRS and good geometry
  if (!grepl("units=m", hybas@proj4string@projargs)){
    hybas <- sp::spTransform(hybas, albers)
  }

  if (!grepl("units=m", station@proj4string@projargs)){
    station <- sp::spTransform(station, albers)
  }

  if (!grepl("units=m", basin@proj4string@projargs)){
    basin <- sp::spTransform(basin, albers)
  }

  P <- rgeos::gBuffer(station, width = 50)

  # Ensure basin is contained within hybas
  basin <- rgeos::gIntersection(basin, hybas)

  # split up 'gaps' between basin and containing hydrobasin
  downstream <- sp::disaggregate(rgeos::gDifference(hybas, basin))
  if (is.null(downstream)){
    warning('Station does not intersect any Hydrobasin - delineation gaps.
            Skipping gap-filling')
    return(basin)
  }

  # Find out which gap contains the station (this one doesn't get filled)
  int <- rgeos::gIntersects(downstream, rgeos::gBuffer(station, width = 100), byid = T)
  if (is.null(nrow(int)) | all(!int)){
    warning('Station does not intersect basin')
    return(basin)
  }

  # Create new geometry with no gaps
  if (additive){
    ## additive approach (fill in gaps)
    final <- rgeos::gUnion(basin, downstream[!int])
  }else{
    ## subtractive approach (subtract whatever is below the station)
    final <- rgeos::gDifference(hybas, downstream[int])
  }

  # try to clean geometries
  final <- outerRing(final)
  final <- invisible(cleangeo::clgeo_Clean(final, strategy = cl_strategy))
  final <- rgeos::gBuffer(final, width = 0)

  #return object similar to original basin object
  if ("data" %in% slotNames(basin)){
    final <- sp::spChFIDs(final, row.names(basin@data))
    final <- sp::SpatialPolygonsDataFrame(final, basin@data)
  }

  if (bprj@projargs != basin@proj4string@projargs){
    final <- sp::spTransform(final, bprj)
  }

  return(final)
}
