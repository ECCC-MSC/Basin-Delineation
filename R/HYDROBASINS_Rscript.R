#' Read Hydrobasins
#'
#' @title Load a HydroBASINS shapefile into R
#' @description  For a directory containing HydroBASIN shapefiles (or for a specified ), load the file
#'    associated with the specified Pfaffstetter level and return an attribute table.
#'    NOTE: If a directory is provided, file naming convention must match original HYBAS (e.g. hybas_na_lev02_v1c.shp)
#' @param directory  the directory containing HydroBASIN shapefiles.  Searches recursively.
#' @param level  The desired Pfaffstetter level
#' @param area  which region to load. can be 'na' (north america), or 'ar' (arctic)
#' @return A list containing a spatial object and an attribute table
#' @examples
#'  LoadHybas("C:\\Data\\HydroBASINS", 'ar', 3)
#'  @export
LoadHybas <- function(path, area='na', level=7){
  if(!file.exists(path)){
  print("File does not exist")
    return()
    }
  if (file.info(path)$isdir == TRUE){
  ## If directory provided, search directory and find desired file (takes first if duplicates)
  match.name <- paste("hybas_", area, sprintf("_lev%02d", level), sep='')
  match.name.reg <- paste(match.name, ".*shp$", sep='')
  files <-  list.files(path, recursive=T, pattern=match.name.reg, ignore.case = T, full.names = T)
  file.to.load <- files[1]

  }else{ #otherwise load file directly
    file.to.load <- path
  }
  ## Read in Shapefile
  print(paste("loading", file.to.load))
  shpdata <- rgdal::readOGR(file.to.load)

  return(shpdata)
}

#' Find Nearest Upstream Sub-Basins
#'
#'@description  For a sub-basin with a given HYBAS.ID, find all sub-basins that are directly upstream (adjacent)
#'uses tables as given from LoadHybas
#'@param HYBAS: a HYBAS spatial polygon object
#'@param HYBAS.ID: The HYBAS.ID of the sub-basin from which to go upstream
#'@param ignore.endorheic Whether or not to ignore 'virtual' inputs from inland basins. defaults to FALSE.
#'@examples
#' HB <- LoadHybas("C:\\Data\\HydroBASINS", 'ar', 3)
#' FindUpstreamSubBasins(HB, "")
#'@export
FindUpstreamSubBasins <- function(HYBAS, HYBAS.ID, ignore.endorheic=F){
  upstream.ab <- HYBAS@data[HYBAS@data$NEXT_DOWN == HYBAS.ID,c("HYBAS_ID", "ENDO")] # returns as char
  if (ignore.endorheic){
    upstream.ab <- upstream.ab[upstream.ab$ENDO != 2,]
  }
  return(as.character(upstream.ab$HYBAS_ID))
}

#' Find All Upstream Sub-Basins
#'
#'@description  For a sub-basin with a given HYBAS.ID, find all sub-basins that are upstream at any distance
#'uses tables as given from LoadHybas
#'@param HYBAS: a HYBAS spatial polygon object
#'@param HYBAS.ID: The HYBAS.ID of the sub-basin from which to go upstream
#'@param ignore.endorheic Whether or not to ignore 'virtual' inputs from inland basins. defaults to FALSE.
#'@export
FindAllUpstreamSubBasins <- function(HYBAS, HYBAS.ID, ignore.endorheic=F, split=F){
  # make containers and stuff
  HYBAS.ID <- as.character(HYBAS.ID)
  HYBAS.ID.master <- list()

  #Get the possible upstream 'branches'
  direct.upstream <- FindUpstreamSubBasins(HYBAS = HYBAS, HYBAS.ID = HYBAS.ID, ignore.endorheic = ignore.endorheic)

  # for each branch iterate upstream until only returning empty results
  for (i in direct.upstream){
    run <- T
    HYBAS.ID.list <- i
    sub.basins <- i # this is the object that gets passed to FindUpstreamSubBasins in each iteration
    while (run){
      result.i <- unlist(lapply(sub.basins, FindUpstreamSubBasins, HYBAS=HYBAS, ignore.endorheic=ignore.endorheic))
      if (length(result.i)==0) run <- F  # Stopping criterion
      HYBAS.ID.list <- c(HYBAS.ID.list, result.i)
      sub.basins <- result.i
    }
    HYBAS.ID.master[[i]] <- HYBAS.ID.list
  }
  if (!split) HYBAS.ID.master <- as.character(unlist(HYBAS.ID.master))

 return(HYBAS.ID.master)
}

#' Find Upstream Stations
#'
#'@description  Finds all HYDAT monitoring stations upstream of given HYBAS sub-basin
#'@param HYBAS: a HYBAS spatial polygon object
#'@param HYBAS.ID: The HYBAS.ID of the sub-basin from which to go upstream
#'@return A character vector of station IDs
#'@keywords internal
#'@export
FindUpstreamStations <- function(HYBAS, HYBAS.ID, con, all=T, ...){
  basin.ID <- FindAllUpstreamSubBasins(HYBAS, HYBAS.ID)
  basins <- HYBAS[HYBAS$HYBAS_ID %in% basin.ID,]
  stations <- SpatialHydat(con)
  upstream.stations <- raster::intersect(stations, basins)
  return(upstream.stations)
}


plot.subset.HYBAS <- function(HYBAS, HYBAS.ID, text=T){
  ## Plots one or more HYBAS polygons
  subset <- HYBAS[HYBAS$HYBAS_ID %in% HYBAS.ID,]
  plot(subset)
  if (text)invisible(text(getSpPPolygonsLabptSlots(subset), labels=as.character(subset$HYBAS_ID), cex=0.7))
}

#' Upstream Basin
#'
#' @description Return the upstream basin polygon for a station
#' @param x either a hybas ID or a HYDAT station ID
#' @param include.first whether or not to include the enclosing hybas polygon (overestimation)
#' @param ignore.endorheic logical, whether or not to include 'virtual connections' in the HydroSHEDS topology
#' @return a spatial polygon data frame of the upstream basin
upstream.basin <- function(station, HYBAS, include.first=T, ignore.endorheic=F){
  station <-NA
  st.ar.gr <- NA
  st.ar.ef <- NA
  if (is.character(x) | is.factor(x)){
    if (grepl('\\d{10}', x)) HYBAS.ID <- as.character(x)
    if (grepl('\\d{2}[A-Za-z]{2}\\d{3}', x)) print("HYDAT needs to be a spatial object")
  }else{
    HYBAS.ID <- as.character(spatialEco::point.in.poly(x,HYBAS)$HYBAS_ID)  # Find ID of containing hydrobasins polygon
    station <-  x[, tolower(names(x)) %in% c("station", "station_number")]
    st.ar.gr <- x[, tolower(names(x)) %in% c("drainage_area_gross", "drainagearea", "drainagear")]
    st.ar.ef <- x[, tolower(names(x)) %in% c("drainage_area_effect")]
    if (length(HYBAS.ID) == 0){ # point outside of HYBAS polygons
      print("Point is outside of HYBAS")  #TODO: This needs a better solution
      out <- NULL
      return(out)
    }
  }
  ## For a given HYBAS ID, get all the upstream basins and their areas
  upstreamIDs <- FindAllUpstreamSubBasins(HYBAS, HYBAS.ID, ignore.endorheic=ignore.endorheic)
  if (include.first){ upstreamIDs <- c(upstreamIDs, HYBAS.ID)}
  subset <- HYBAS[HYBAS$HYBAS_ID %in% upstreamIDs,]

  # Transform to get area
  subset.t <-  rgeos::gUnaryUnion(sp::spTransform(subset, CRSobj = CRS.AlbersEqualAreaConic))
  area <- rgeos::gArea(subset.t, byid = F)/1e6

  # Build a filled area -  TODO: this is possibly outdated.  can probably delete.
  # n <- which.max(sapply(subset.t@polygons[[1]]@Polygons, slot, name='area'))
  # ring = SpatialPolygons(list(Polygons(list(
  #  subset.t@polygons[[1]]@Polygons[[n]]),
  # ID=1)))
  # ring@proj4string <- subset.t@proj4string
  # area.ring <- rgeos::gArea(ring, byid = F)/1e6

  # merge polygon and store area measurements
  out <- rgeos::gUnaryUnion(subset)
  row.names(out) <- as.character(1:length(out))
  data <- as.data.frame(list(
    HYBAS=HYBAS.ID,
    station=station,
    station_eff=st.ar.ef,
    station_gro=st.ar.gr,
    area=area,
    filled_area=area.ring,
    HYB_UP_AR=HYBAS@data[HYBAS$HYBAS_ID==HYBAS.ID,"UP_AREA"]))
  out <- SpatialPolygonsDataFrame(out, data)
  return(out)
  #qur <- list("sp.polygons", querypoly, fill='red')
  #spplot(sp.layout=list(ups, qur))
}

#' Upstream Basin 2
#'
#' @description Return the upstream basin polygon for a station
#' @param station A SpatialPointsDataFrame corresponding to the station of intrest
#' @param ignore.endorheic logical, whether or not to include 'virtual connections' in the HydroSHEDS topology
#' @return a spatial polygon data frame of the upstream basin in the CRS of the hybas polygon (should be WGS84)
UpstreamHYBAS <- function(station, HYBAS, ignore.endorheic=F, fast=T){

  if (!station@proj4string@projargs == HYBAS@proj4string@projargs){  # Ensure both are in HYBAS CRS (likely WGS84)
    station <- sp::spTransform(station, CRSobj = CRS(HYBAS@proj4string@projargs))
  }

  if (fast){ # subset basins first
    HYBAS <- HYBAS[substr(HYBAS$PFAF_ID,1,3) %in% NearestHYBAS(station@data$station_number), ]
  }

  HYBAS.ID <- as.character(spatialEco::point.in.poly(station, HYBAS)$HYBAS_ID)  # Find ID of containing hydrobasins polygon
    if (length(HYBAS.ID) == 0){ # point outside of HYBAS polygons
      print("Point is outside of HYBAS")  #TODO: This needs a better solution
      return(NULL)
    }else{
      print(sprintf("Station is within HYBAS polygon %s", HYBAS.ID))
    }

  ## For a given HYBAS ID, get all the upstream basins and their areas
  upstreamIDs <- FindAllUpstreamSubBasins(HYBAS, HYBAS.ID, ignore.endorheic=ignore.endorheic, split = T)
  out <- HYBAS[match(c(HYBAS.ID ,unlist(upstreamIDs)), HYBAS@data$HYBAS_ID, nomatch=0),]
  id <- list()
  if (length(upstreamIDs)!=0){
    names(upstreamIDs) <- LETTERS[2:(length(upstreamIDs)+1)]
    id <- lapply(names(upstreamIDs), function(x) upstreamIDs[[x]] <- rep(x, length(upstreamIDs[[x]])))
  }
  out@data$ID <- c("A",unlist(id))
 # out <- rgeos::gUnaryUnion(out, id=c("A",unlist(id)))
  return(out)
}


#' Find nearest river segment
#'
#' @description Find the nearest line segment (river) to a point (station).  Spatial objects should be in an
#' appropriate coordinate system (i.e. not lat-long) otherwise they will be projected and this adds to processing time.
#' If CRS does not match between the stations and the rivers they will both be projected into Canada Albers Equal
#' Area conic.
#' @param spatialPoint A SpatialPointsDataFrame of the desired station
#' @param riverLines A SpatialLinesDataFrame of the HydroSHEDS river network
#' @return dataframe of nearest lines (attribute table of matches)
FindNearestRiverSegment <- function(spatialPoint, riverLines, HYBAS){
  if (!missing(HYBAS)){
    # (optional) clip rivers using hybas to speed up transformation
    riverLines <- gIntersection(HYBAS, rivers, byid = TRUE, drop_lower_td = TRUE)
    }

  # reproject rivers and point into UTM if necessary
  if (spatialPoint@proj4string@projargs != riverLines@proj4string@projargs){
    utm.zone <-  WhichZone(spatialPoint@data$longitude)
    CRS.string <- paste("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs", sep='')#+proj=utm +zone=",utm.zone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs
    spatialPoint <- sp::spTransform(spatialPoint, CRSobj = CRS(CRS.string))
    HYBAS <- sp::spTransform(HYBAS, CRSobj = CRS(CRS.string))
  }


  # find nearest river segment (?)
  result <- data.frame(matrix(numeric(0), ncol=2))
  distance <- numeric()
  names(result) <- names(riverLines@data)
  for (i in seq(1:nrow(spatialPoint@data))){
    delta <- rgeos::gDistance(spatialPoint[i,], riverLines, byid=TRUE)
    smallest <- which.min(delta)
    x <- riverLines[smallest,]
    result <- rbind(result,x@data)
    distances <- c(distances, delta[smallest])
    snappedcoords <- t(as.matrix(maptools::nearestPointOnLine(
      x@lines[[1]]@Lines[[1]]@coords, spatialPoint[i,]@coords)))
  }

  #return coordinate of point
 # maptools::nearestPointOnLine(x@lines[[1]]@Lines[[1]]@coords, spatialPoint[i,]@coords )
  # return ID of nearest river segment(s)
  return(cbind(result, distances, snappedcoords))
}


#' Process Station
#'
#' @description Gives information about a station (upstream area, upstream stations)
#' @param station Either a SpatialPointsDataFrame or the path to a Hydat-exported table with desired points
#' or a character vector of station names (in this case a connection must be specified)
#' @param con (optional) A open database connection.  Used if station is a character vector of station names
#' @param DEM either
#' @param outdir (optional) A directory to which the output shapefile and report will be saved.  If not specified,
#' the function will return an R object without saving data to disk.
#' @param pourpoints spatial points data frame with
#' @param NCA (optional) Either a character string path or a SpatialPolygonsDataFrame.  This is a polygon
#' of the non-contributing areas.  If provided, calculates effective drainage area
#' @return A SpatialPolygonDataFrame
StationPolygon <- function(station, con, outdir, HYBAS, DEM.path, DEM.source, saga.env, NCA, interactive=F,pourpoints,
                           projected.CRS="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs", ...){

  # Figure out what "station" is and read it in
  if (class(station) == "SpatialPointsDataFrame"){
   print("detected a spatial object for: 'station'")
  }else if (!missing(con) & class(station)=="character"){
    print("reading 'station' from database connection")
    station <- SpatialHydat(con, station)
  }else if (file.exists(station)){
    print("reading 'station' from file")
    station <- SpatialECDE(station)
  }else{
    print("could not read 'station' information")
  }
  print(station$station_number)
  if (!missing(pourpoints)){
    station <- getPourPoint(station, pourpoints)
  }
  #=============
  # Calculations
  #=============

  # Get upslope area from HYBAS
  HYB.Poly.individual <- UpstreamHYBAS(station, HYBAS)
  HYB.Poly <- rgeos::gUnaryUnion(HYB.Poly.individual, id=HYB.Poly.individual$ID) # combined area
  no.upstream.hyb <- (length(HYB.Poly.individual) == 1)

  # Get first two upstream tiles (should be the first one for each ID)
  first.up <- unlist(by(data = HYB.Poly.individual@data, INDICES = HYB.Poly.individual@data$ID,
                     FUN = function(x) as.character(x[which.max(as.character(x$UP_AREA)),c("HYBAS_ID")])))
  first.up.HYB <- HYB.Poly.individual[HYB.Poly.individual$HYBAS_ID %in% first.up & HYB.Poly.individual$ENDO != 2,]
  first.up.HYB.endo <- HYB.Poly.individual[HYB.Poly.individual$HYBAS_ID %in% first.up,]


  # find any endorheic basins
  endo.HYB <- HYB.Poly.individual@data[HYB.Poly.individual$HYBAS_ID %in% first.up & HYB.Poly.individual$ENDO == 2,"ID"]
  ID.list <- sapply(slot(HYB.Poly, "polygons"), function(x) slot(x, "ID"))
  endo <- ID.list %in% endo.HYB

  disjoint <- rep(F, length(endo))
  not.unlikely <- rep(F, length(endo))

  # find necessary coverage for DEM

  DEM.path <- OverlayNTS(gUnaryUnion(first.up.HYB),
                         NTS50k, saga.env$workspace, saga.env$workspace,
                         tilename="NTS_SNRC")

  # Get upslope area from DEM
  DEMresult <- UpslopeDEM(station, DEM.path = DEM.path, saga.env = saga.env,
                         outdir = saga.env$workspace, projected.CRS=projected.CRS, ...)
  if (DEMresult$nodata==T){ # if the point has no data
    DEM.poly.p <- sp::spTransform(HYB.Poly, CRS = sp::CRS(projected.CRS)) #use different CRS?
  }
  DEM.Poly <- DEMresult$final.name
  pointbuffer <- DEMresult$pointbuffer
  iterate.to <- DEMresult$iterate.to
  iterate.thres <- DEMresult$iterate.thres

  DEM.Poly.p <- invisible(rgdal::readOGR(DEM.Poly))
  DEM.Poly.p <- DEM.Poly.p[which.max(area(DEM.Poly.p)),] # take the big one.  others will likely have bad geometries


  # Clip 'nose' of basin
  HP.p <- sp::spTransform(HYB.Poly, CRSobj = sp::CRS(DEM.Poly.p@proj4string@projargs))
  local.drainage.p <- gIntersection(HP.p[1,], DEM.Poly.p, byid = F, drop_lower_td = TRUE) # DEM area within containing hybas


  #=====================
  #  Logical Steps
  #=====================

  if (!no.upstream.hyb){
    # if DEM is disjoint with upslope area, do not include it.
    disjoint <- gDisjoint(HP.p, local.drainage.p, byid = T)

    # if calculated total upslope area is less than half of an upslope HYBAS region, remove that HYBAS region
    not.unlikely <- area(DEM.Poly.p) > 0.5*area(sp::spTransform(first.up.HYB.endo, CRS(DEM.Poly.p@proj4string@projargs)))
  }

  keep <- ID.list[(!disjoint | endo) & (not.unlikely | endo)]
  keep <- keep[-which(keep=="A")]
  HP.out <- HP.p[keep,]

  # merge desired parts of HYBAS polygons with local drainage
  if (pointbuffer >= iterate.to){
    output <- gUnaryUnion(HP.p)
  }else if (length(HP.out) > 0 & !no.upstream.hyb){
    HP.p <- gUnaryUnion(HP.out)
    output <- gUnion(local.drainage.p, HP.out)
  }else{
    output <- gUnaryUnion(local.drainage.p)
  }

  dr_ar <- area(output) *1e-6

  # clean geometry with a zero-width buffer
  output <- gBuffer(output, byid=TRUE, width=0)

  # Calculate non-contributing areas
  if (!missing(NCA)){
    NCA <- InterpretShapefile(NCA)
    if (NCA@proj4string@projargs != projected.CRS){
    NCA <- sp::spTransform(NCA, sp::CRS(projected.CRS))
    }
     clipped <- gDifference(output, NCA)
     drainage.effective <- round(area(clipped)) * 1e-6
    print("take out non-contributing area") # take out non-contributing area
  }else{
    drainage.effective <- NA
  }

  # create attribute table
  data <- as.data.frame(list(  # make output table
    stn_number=station@data$station_number,
    drainar_gr=dr_ar,
    drainar_ef=drainage.effective,  # add notes?
    stn_dr_gr=station@data$drainage_area_gross,
    stn_dr_ef=station@data$drainage_area_effect,
    trgtRdius=pointbuffer,
    trgtRdMax=iterate.to,
    iterThres=iterate.thres,
    snapped=DEMresult$snapped,
    snap.dist=DEMresult$snap.dist
  ))

   output <-  sp::spTransform(SpatialPolygonsDataFrame(output, data), # output shape file with same CRS as station
                              CRSobj = sp::CRS(station@proj4string@projargs))

  # write output shapefile (rgdal::writeOGR)
  name <- paste(station$station_number, "_basin.shp")
  rgdal::writeOGR(output, dsn=file.path(outdir, name), layer=name, driver="ESRI Shapefile")

  # write output table (optional?)
  tablename <- paste(station$station_number, "_basin.csv")
  write.csv(output@data, file.path(outdir, tablename), row.names=F)

  #clean up (optional?)
  file.remove(DEM.Poly)

  return(list(hyb = HP.p, dem = local.drainage.p, out=output))
}

#' Upstream Area
#'
#' @description Finds basin area from raster flow direction and accumulation maps
#' @param point spatial point or polygon
#' @param flow.direction
#' @param flow.accu flow accumulation raster
UpstreamBasin <- function(point, flow.direction, flow.accu){
return(NULL)
}

# get pour point
getPourPoint <- function(point, pourpoint){
  if (point@data$station_number %in% pourpoint@data$station_number) {
    pp <- pourpoint[pourpoint@data$station_number == point@data$station_number, "station_number"]
    pp@data <- merge(pp@data, point@data)
    print('Using pour point coordinate')
    return(pp)
  }else{
    print('no matching pour point found!')
    return(point)
  }
}
