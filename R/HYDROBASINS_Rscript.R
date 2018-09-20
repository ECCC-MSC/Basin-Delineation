#===============================================================================
#' @title Load a HydroBASINS shapefile into R
#'
#'
#' @description  For a directory containing HydroBASIN shapefiles , load the file
#' associated with the specified Pfaffstetter level and return an attribute table.
#'NOTE: If a directory is provided, file naming convention must match original
#'HYBAS (e.g. hybas_na_lev02_v1c.shp)
#'
#' @param directory  the directory containing HydroBASIN shapefiles.
#' Searches recursively.
#'
#' @param level  The desired Pfaffstetter level
#'
#' @param area  which region to load. can be 'na' (north america), or 'ar'
#'  (arctic)
#'
#' @return A list containing a spatial object and an attribute table
#'
#' @examples
#'  LoadHybas("~/HydroBASINS", 'ar', 3)
#'
#' @export
#===============================================================================
LoadHybas <- function(path, area = "na", level = 7){
  if (!file.exists(path)){
  print("Directory does not exist")
    return()
  }

  if (file.info(path)$isdir == TRUE){

  ## If directory provided, search directory and find desired file (takes first if duplicates)
  match_name <- paste("hybas_", area, sprintf("_lev%02d", level), sep = "")
  match_name_reg <- paste(match_name, ".*shp$", sep = "")
  files <-  list.files(path, recursive = T, pattern = match_name_reg,
                       ignore.case = T, full.names = T)

  file_to_load <- files[1]

  }else{
    # if not a directory, load file directly
    file_to_load <- path
  }

  ## Read in Shapefile
  print(paste("loading", file_to_load))
  shpdata <- rgdal::readOGR(file_to_load, stringsAsFactors = FALSE)

  return(shpdata)
}

#===============================================================================
#' @title Find Downstream Basins
#'
#'
#' @description Uses the HydroBasins "NEXT_DOWN" field to identify all downstream
#' polygons from a target polygon. It is possible to stop downstream navigation
#' at endorheic sinks (terminal polygons within inland basins), or to continue
#' through 'virtual' connections (see HydroBASINS technical documentation for
#'  more details).
#'
#' @param HYBAS.ID character vector of one or more HydroBASINS ID's
#'
#' @param HYBAS a HYBAS spatial polygon object or dataframe from such an object
#'
#' @param n_iter integer, how many iterations to perform. To remove limit,
#' set to -1 (default)
#'
#' @param ignore.endo logical, whether or not to include 'virtual' connections
#' from endorheic sinks (see HydroBASINS technical documentation for more details)
#'
#' @return a character vector of HydroBASIN ID's corresponding to downstream
#' basins
#===============================================================================
FindDownstreamBasins <- function(HYBAS, HYBAS.ID, n_iter = -1, ignore.endo = T){
  HYBAS.NEW <- HYBAS.ID
  i <- 0
  while (length(HYBAS.NEW) > 0 & (n_iter == -1 | i < n_iter)){

    if (ignore.endo){
      HYBAS.DOWN <- HYBAS@data[HYBAS$HYBAS_ID %in% HYBAS.NEW & HYBAS$ENDO != 2,
                               "NEXT_DOWN"]
    }else{
      HYBAS.DOWN <- HYBAS@data[HYBAS$HYBAS_ID %in% HYBAS.NEW, "NEXT_DOWN"]
    }

    HYBAS.DOWN <- HYBAS@data[HYBAS$HYBAS_ID %in% HYBAS.DOWN, "HYBAS_ID"]
    HYBAS.NEW <- HYBAS.DOWN[!(HYBAS.DOWN %in% HYBAS.ID)]
    HYBAS.ID <- c(HYBAS.ID, HYBAS.NEW)

    i <- i + 1
  }

  return(HYBAS.ID[-1])
}


#===============================================================================
#' @title Find Nearest Upstream Sub-Basins
#'
#'
#' @description  For a sub-basin with a given HYBAS.ID, find all sub-basins
#' that are directly upstream (adjacent)
#' uses tables as given from LoadHybas
#'
#' @param HYBAS: a HYBAS spatial polygon object or dataframe from such an object
#'
#' @param HYBAS.ID: The HYBAS.ID of the sub-basin from which to go upstream
#'
#' @param ignore.endorheic Whether or not to ignore 'virtual' inputs from
#' inland basins. defaults to FALSE.
#'
#' @examples
#' HB <- LoadHybas("C:\\Data\\HydroBASINS", 'ar', 3)
#' FindUpstreamSubBasins(HB, "")
#'
#' @export
#===============================================================================
FindUpstreamSubBasins <- function(HYBAS, HYBAS.ID, ignore.endorheic = F){
  if (!class(HYBAS) == "data.frame"){
    HYBAS <- HYBAS@data
  }

  upstream.ab <- HYBAS[HYBAS$NEXT_DOWN == HYBAS.ID, c("HYBAS_ID", "ENDO")] # returns as char
  if (ignore.endorheic){
    upstream.ab <- upstream.ab[upstream.ab$ENDO != 2, ]
  }
  return(as.character(upstream.ab$HYBAS_ID))
}

#===============================================================================
#' @title Find All Upstream Sub-Basins
#'
#'
#' @description  For a sub-basin with a given HYBAS.ID, find all sub-basins that
#' are upstream at any distance
#' uses tables as given from LoadHybas
#'
#' @param HYBAS: a HYBAS spatial polygon object or dataframe from such an object
#'
#' @param HYBAS.ID: The HYBAS.ID of the sub-basin from which to go upstream
#'
#' @param ignore.endorheic Whether or not to ignore 'virtual' inputs from inland
#'  basins. defaults to FALSE.
#'
#' @export
#===============================================================================
FindAllUpstreamSubBasins <- function(HYBAS, HYBAS.ID, ignore.endorheic = F,
                                     split = F){
  # Simplify HYBAS to a dataframe
  if (!class(HYBAS) == "data.frame"){
    HYBAS <- HYBAS@data
  }

  # make containers
  HYBAS.ID <- as.character(HYBAS.ID)
  HYBAS.ID.master <- list()

  #Get the possible upstream 'branches'
  direct.upstream <- FindUpstreamSubBasins(HYBAS = HYBAS, HYBAS.ID = HYBAS.ID,
                                           ignore.endorheic = ignore.endorheic)

  # for each branch iterate upstream until only returning empty results
  for (i in direct.upstream){
    run <- T
    HYBAS.ID.list <- i
    sub.basins <- i # this is the object that gets passed to FindUpstreamSubBasins in each iteration
    while (run){
      result.i <- unlist(lapply(sub.basins, FindUpstreamSubBasins,
                                HYBAS = HYBAS, ignore.endorheic = ignore.endorheic))

      if (length(result.i) == 0){ run <- F } # Stopping criterion
      HYBAS.ID.list <- c(HYBAS.ID.list, result.i)
      sub.basins <- result.i
    }
    HYBAS.ID.master[[i]] <- HYBAS.ID.list
  }

  if (!split){ HYBAS.ID.master <- as.character(unlist(HYBAS.ID.master)) }

 return(HYBAS.ID.master)
}


#===============================================================================
#' Create HydroBASINS catchment basin
#'
#'
#' @description Return the upstream basin polygon for a station
#'
#' @param station A SpatialPointsDataFrame corresponding to the station of
#' intrest
#'
#' @param ignore.endorheic logical, whether or not to include 'virtual
#' connections' in the HydroSHEDS topology
#'
#' @param outputfile (optional) if provided, saves the upstream hydrobasins to
#' a shapefile. File name must
#' end with '.shp'
#'
#' @return a spatial polygon data frame of the upstream basin in the CRS of the
#'  hybas polygon (should be WGS84)
#'
#' @export
#===============================================================================
UpstreamHYBAS <- function(station, HYBAS, ignore.endorheic = F, fast = T,
                          outputfile){

  if (!station@proj4string@projargs == HYBAS@proj4string@projargs){
    # Ensure both are in HYBAS CRS (WGS84)
    station <- sp::spTransform(station, CRSobj = sp::CRS(HYBAS@proj4string@projargs))
  }

  if (fast){ # subset basins first to improve speed
    HYBAS <- HYBAS[substr(HYBAS$PFAF_ID, 1, 3) %in%
                     NearestHYBAS(station@data$station_number), ]
  }

  HYBAS.ID <- as.character(spatialEco::point.in.poly(station, HYBAS)$HYBAS_ID)  # Find ID of containing hydrobasins polygon
    if (length(HYBAS.ID) == 0){ # point outside of HYBAS polygons
      print("Point is outside of HYBAS")  #TODO: This needs a better solution
      return(NULL)
    }else{
      print(sprintf("Station is within HYBAS polygon %s", HYBAS.ID))
    }

  ## For a given HYBAS ID, get all the upstream basins and their areas
  upstreamIDs <- FindAllUpstreamSubBasins(HYBAS, HYBAS.ID,
                                          ignore.endorheic = ignore.endorheic,
                                          split = T)

  out <- HYBAS[match(c(HYBAS.ID, unlist(upstreamIDs)), HYBAS@data$HYBAS_ID,
                     nomatch = 0), ]

  id <- list()
  if (length(upstreamIDs) != 0){
    names(upstreamIDs) <- LETTERS[2:(length(upstreamIDs) + 1)]
    id <- lapply(names(upstreamIDs),
               function(x) upstreamIDs[[x]] <- rep(x, length(upstreamIDs[[x]])))
  }
  out@data$ID <- c("A", unlist(id))
  if (!missing(outputfile)){
    rgdal::writeOGR(out, dsn = outputfile, layer = basename(outputfile),
                    driver = "ESRI Shapefile")
  }

  return(out)
}

#===============================================================================
#' @title Find nearest river segment
#'
#' @description Find the nearest line segment (river) to a point (station).
#' Spatial objects should be in an appropriate coordinate system (i.e. not
#' lat-long) otherwise they will be projected and this adds to processing time.
#' If CRS does not match between the stations and the rivers they will both be
#' projected into Canada Albers Equal Area conic.
#'
#' @param spatialPoint A SpatialPointsDataFrame of the desired station
#'
#' @param riverLines A SpatialLinesDataFrame of the HydroSHEDS river network
#'
#' @return dataframe of nearest lines (attribute table of matches)
#===============================================================================
FindNearestRiverSegment <- function(spatialPoint, riverLines, HYBAS){
  if (!missing(HYBAS)){
    # (optional) clip rivers using hybas to speed up transformation
    riverLines <- rgeos::gIntersection(HYBAS, rivers, byid = TRUE,
                                       drop_lower_td = TRUE)
    }

  # reproject rivers and point into UTM if necessary
  if (spatialPoint@proj4string@projargs != riverLines@proj4string@projargs){
    CRS.string <- paste("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96
                        +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m
                        +no_defs", sep = "")
    spatialPoint <- sp::spTransform(spatialPoint, CRSobj = sp::CRS(CRS.string))
    HYBAS <- sp::spTransform(HYBAS, CRSobj = sp::CRS(CRS.string))
  }


  # find nearest river segment (?)
  result <- data.frame(matrix(numeric(0), ncol = 2))
  distances <- numeric()
  names(result) <- names(riverLines@data)

  for (i in seq(1:nrow(spatialPoint@data))){
    delta <- rgeos::gDistance(spatialPoint[i, ], riverLines, byid = TRUE)
    smallest <- which.min(delta)
    x <- riverLines[smallest, ]
    result <- rbind(result, x@data)
    distances <- c(distances, delta[smallest])
    snappedcoords <- t(as.matrix(maptools::nearestPointOnLine(
      x@lines[[1]]@Lines[[1]]@coords, spatialPoint[i, ]@coords)))
  }

  #return coordinate of point
 # maptools::nearestPointOnLine(x@lines[[1]]@Lines[[1]]@coords, spatialPoint[i,]@coords )
  # return ID of nearest river segment(s)
  return(cbind(result, distances, snappedcoords))
}

#===============================================================================

#===============================================================================
ReadStationInformation <- function(station, con){
  # Figure out what "station" is and read it in
  if (class(station) == "SpatialPointsDataFrame"){
    print("detected a spatial object for: 'station'")
  }else if (!missing(con) & class(station) == "character"){
    print("reading 'station' from database connection")
    station <- SpatialHydat(con, station)
  }else if (file.exists(station)){
    print("reading 'station' from file")
    station <- SpatialECDE(station)
  }else{
    print("could not read 'station' information")
  }
  print(station$station_number)
  if (nrow(station) == 0){stop("'station' object contains no points.")}
  return(station)
}

#===============================================================================
#' @title Calculate non-contributing area for
#'
#===============================================================================
CalculateNCA <- function(NCA, basin, projected.CRS, verbose=T){
  if (verbose){cat("calculating non-contributing area ... ")}
  NCA <- InterpretShapefile(NCA)
  if (NCA@proj4string@projargs != projected.CRS){
    NCA <- sp::spTransform(NCA, sp::CRS(projected.CRS))
  }
  basin <- InterpretShapefile(basin)
  if (basin@proj4string@projargs != projected.CRS){
    basin <- sp::spTransform(basin, sp::CRS(projected.CRS))
  }
  clipped <- rgeos::gDifference(basin, NCA)
  drainage_effective <- round(raster::area(clipped)) * 1e-6

  if (verbose){cat(sprintf(" %.1f km2 \n", drainage_effective))}
  return(drainage_effective)
}


#===============================================================================
#' @title Move station to pourpoint
#'
#' @description changes the location of a station (point) to that of its
#' pourpoint using the 'station_number' column as a shared ID.
#'
#' @param point a SpatialPointsDataFrame representing one station. Must have a
#' column named 'station_number'
#'
#' @param pourpoint a SpatialPointsDataFrame of pourpoints. Must have a column
#' named 'station_number'
#'
#' @return A SpatialPointsDataFrame with data identical to the original point,
#' and the new geometry corresponding to the pour point. If no matching pour
#' point is found, the original station point is returned
#===============================================================================
getPourPoint <- function(point, pourpoint){
  if (point@data$station_number %in% pourpoint@data$station_number) {
    pp <- pourpoint[pourpoint@data$station_number ==
                      point@data$station_number, "station_number"]
    pp@data <- merge(pp@data, point@data)
    print("Using pour point coordinate")
    return(pp)
  }else{
    print("No matching pour point found!")
    return(point)
  }
}
