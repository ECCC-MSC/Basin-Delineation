
#===============================================================================
#' @title Create containing HydroBASINS
#'
#' @description Uses a classification code to build the maximum allowable
#'  drainage basin for a station using HydroBASINS polygons.
#'
#' @param df data.frame with exactly 1 row containing station names and the
#'  HYBAS_ID of the HydroBASINS polygon that contains each station.
#'
#' @param HYBAS SpatialPolygonsDataFrame of hydrobasins polygons
#'
#' @param code character, classification code describing location of station
#' within cell and on stream
#'
#' @param outdir character, directory in which to save output shapefile
#'
#' @param hyb.id character, name of column containing ID of HydroBASINS sub-basin
#'  that contains the station
#'
#' @param stn.id character, name of column containing the station number
#'
#' @details uses the following classification system
#' P: containing sub-basin has 2 input sub-basins
#'   Pcb Waterbody continues upstream and station on main branch below confluence
#'   Pct Waterbody continues upstream and station on tributary branch above confluence
#'   Pcm Waterbody continues upstream and station on main branch above confluence
#'   Ps Waterbody terminates upstream within sub-basin
#'
#' I: containing sub-basin has 1 input sub-basins
#'   Ic Waterbody continues upstream
#'   Is Waterbody terminates upstream within sub-basin
#'
#' T: containing sub-basin has 0 input sub-basins
#'
#' The output shapefile contains the containing hydrobasin.
#' @export
#===============================================================================
HYBASBasinLimits <- function(df, HYBAS, code, outdir,
         hyb.id='HYBAS_ID', stn.id = 'station_number'){
  print(paste(df[,stn.id], code))

  # get upstream HYBAS
  if (code %in% c("Pcb", "Ic")){
    ids <- FindAllUpstreamSubBasins(HYBAS=HYBAS, HYBAS.ID = df[,hyb.id])

  }else if (code %in% c("Pct")){
    next2up <- FindUpstreamSubBasins(HYBAS=HYBAS,
                                     HYBAS.ID = df[,hyb.id],
                                     ignore.endorheic = T)

    next2up <- HYBAS[HYBAS$HYBAS_ID %in% next2up,]
    nextup <- next2up@data[which.min(next2up$UP_AREA), "HYBAS_ID"] # find ID of next-up on smaller branch
    ids <- FindAllUpstreamSubBasins(HYBAS=HYBAS, HYBAS.ID = nextup)
    ids <- c(ids, nextup)

  }else if (code %in% c("Pcm")){
    next2up <- FindUpstreamSubBasins(HYBAS=HYBAS, HYBAS.ID = df[,hyb.id],
                                     ignore.endorheic = T)

    next2up <- HYBAS[HYBAS$HYBAS_ID %in% next2up,]
    nextup <- next2up@data[which.max(next2up$UP_AREA), "HYBAS_ID"] # find ID of next-up on main branch
    ids <- FindAllUpstreamSubBasins(HYBAS=HYBAS, HYBAS.ID = nextup)
    ids <- c(ids, nextup)

  }else if (grepl("[Ts]", code) & !grepl("du", code)){
    ids <- character()

  }else{
    print("not a valid code")
    return()
  }
  ids <- c(df[,hyb.id], ids)

  # save output file
  out <- HYBAS[HYBAS$HYBAS_ID %in% ids,]
  outname <- paste(df[,stn.id], "HYBAS_Basin.shp", sep='_')
  outfile <- file.path(outdir, outname)
  writeOGR(obj = out, dsn = outfile, layer = gsub('\\.shp','', outname),
           driver = "ESRI Shapefile", overwrite_layer = T)
}

#===============================================================================
#' @title Combine Watersheds
#'
#' @description Combine DEM-derived and HydroBASINS-derived watershed
#'
#' @param station station number in format "01XY001"
#'
#' @param shedspoly_folder character, path to folder containing maximum
#'  drainage basin shapefiles derived from a HydroSHEDS
#'
#' @param dempoly_folder character, path to folder containing drainage basin
#' shapefiles derived from DEM
#'
#' @param out_folder character, path to folder containing output files
#'
#' @export
#===============================================================================
CombineWatersheds <- function(station, shedspoly_folder, dempoly_folder,
                              out_folder){
  if (!class(station)=='character'){stop("station should be a character string")}
  # find hydrobasins polygon in folder
  HYB_shp <- list.files(shedspoly_folder, pattern = paste0(station, ".*\\.shp$"),
                        full.names = T)
  if (length(HYB_shp)>1){warning("Multiple hydrobasins polygons found. Using the
                                 first one.")}
  HYB <- rgdal::readOGR(HYB_shp[1])

  # find DEM polygon in folder
  DEM_shp <- list.files(dempoly_folder, pattern = paste0(station, ".*\\.shp$"),
                        full.names = T)
  if (length(DEM_shp)>1){warnings("Multiple DEM polygons found. Using the
                                  first one.")}
  DEM <- rgdal::readOGR(DEM_shp[1])
  DEM_data_csv <- list.files(dempoly_folder, pattern = paste0(station, ".*\\.csv$"),
                             full.names = T)
  DEM_data <- read.csv(DEM_data_csv)

  # reproject DEM polygon
  DEM <- sp::spTransform(DEM, CRSobj = CRS(HYB@proj4string@projargs))


  # Perform set operations
  A = HYB[-which.max(HYB$UP_AREA),]
  C = HYB[which.max(HYB$UP_AREA),]

  output <- rgeos::gIntersection(DEM, HYB)
  if (nrow(A)!=0){output <- rgeos::gUnion(output, A)}

  # remove 'floating' polygons
  output <- sp::disaggregate(output)
  output <- output[which.max(lapply(output@polygons, slot, name='area')),]

  # remove holes (inner rings)
  output <- outerRing(output)

  # add in data to attributes table
  data <- DEM_data[,c('station_number', 'pointbuffer', 'nodata', 'snap.dist',
                      'DEM', 'cell_acc', 'method')]
  data$area <- round(raster::area(output)*1e-6,2)
  rownames(data) <- sapply(output@polygons, slot, name='ID')
  output <- SpatialPolygonsDataFrame(output, data)

  output <- sp::spTransform(output,
              CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

  # write output
  file_out <- file.path(out_folder, paste0(station,"_cntrb_basin.shp"))
  rgdal::writeOGR(obj = output,
           dsn = file_out,
           layer = basename(file_out),
           driver = 'ESRI Shapefile')
  write.csv(output@data, sub("shp$", "csv", file_out), row.names=F, quote=F)
  }

#===============================================================================
#' @title Calculate basin area using a DEM
#'
#' @description Create a drainage basin shapefile using a DEM and a
#' @param point a SpatialPointsDataFrame corresponding to a hydrometric station.
#'  Must be in same coordinates system as DEM and coordinate system
#' must be projected (not lat/long).  Must have longitude and latitude attributes.
#'
#' @param saga.env an rsaga environment object
#'
#' @param outdir Directory to output final upslope shapefile
#'
#' @param pointbuffer numeric, how much should the point be buffered (calculates
#' upslope basin of buffered area)
#'
#' @param DEM.source character, one of: c('CDED', 'NED', 'CDEM, CDSM', 'SHEDS').
#'  Ignored if a DEM file is supplied to DEM.path
#'
#' @param DEM.path One of:
#'  (1) a file path to a directory containing DEM files either in the the format
#'   n%%w0%%_con_grid.sgrd (if DEM.source is 'SHEDS') or a directory to which DEM
#'   tiles will be downloaded (if DEM.source='CDEM', 'CDED' etc.),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a
#'  projected coordinate system and the coordinate system should match that of
#'  the point (e.g. Canada Albers Conformal Conic)
#'
#' @param dem.clip.square How big (in m) a clip should be generated from the
#' original DEM (too big doesn't work and is slower)
#'
#' @param projected.CRS character vector (proj4) specifying "working" CRS in
#' which area and distance calculations are to be done
#'
#' @param method character, one of ('D8', 'DINF'), specifying the hydrological
#' pathing model
#'
#' @param upstream.threshold integer, number of upstream cells required to
#' initiate channel growth. A smaller number will yield more channels and will
#'  snap to smaller streams. Larger numbers will produce only major channel
#'  branches. For a 50k DEM, numbers between 20,000 and 50,000 are appropriate
#'  to reproduce the main channel network without many small branches. Numbers
#'  around 8,000 will produce most small branches, but tend to create many offshots
#'  in wide river sections.
#'
#' @return A list containing diagnostic information and and a character string
#' naming the recently created upstream area polygon
#'
#' @export
#===============================================================================
DEMDrainageBasin <- function(point, DEM.path, DEM.source='NTS', saga.env,
                             outdir, pointbuffer=50, dem.clip.square=50000,
                             projected.CRS,
                             upstr_method="D8", upstream.threshold=35000, ...){

  # set some parameters
  method <- switch(toupper(upstr_method), "D8"=0, "DINF"=1)
  snapped<-F
  pb.init <- pointbuffer

  #sanity checks
  if (missing(projected.CRS)){
    projected.CRS <- GetProj4("AlbersEqualAreaConic")
  }
  in.DEM <- gsub("\\.sdat$", "\\.sgrd", DEM.path)

  # prepare workspace
  oldwd <- getwd()
  setwd(saga.env$workspace)
  outdir <- gsub("[\\/]$","",outdir)

  # get name of station
  station.name <- point@data[1,grepl("^\\d{2}[[:alpha:]]{2}\\d{3}$", point@data)]
  print(station.name)
  final.name <- file.path(outdir, paste(station.name,"_upslope.shp",sep=''))

  #if point not in proper crs, convert it
  if (point@proj4string@projargs != projected.CRS){
    point <- sp::spTransform(point, CRSobj = CRS(projected.CRS))
  }

  # Find missing DEM if necessary / Convert to sgrd if necessary.
  if (!grepl("\\.sgrd$", in.DEM)){  # if not an SGRD file (then we expect a directory or a list of names)

    if (tolower(DEM.source)=='sheds'){
      name <- HydroMosaic(point@data$longitude, point@data$latitude,
                          tol=dem.clip.square)

      if (length(name)==1){
        check.if.exists <- list.files(in.DEM,
                                      pattern=paste(name, ".sgrd", sep=''),
                                      full.names = T, recursive = T)

        if (length(check.if.exists) ==0){  #  if sgrid doesn't exist, create it
          original.DEM <- GetTilePathsHS(name, in.DEM)
          print(sprintf("Converting %s DEM to sgrd...", original.DEM))
          dstfile <- gsub("_con$","_con\\.sdat",original.DEM)
          gdal_warp2SAGA(original.DEM, outputCRS = projected.CRS, dstfile = dstfile)
        }

        in.DEM <- list.files(in.DEM, pattern=paste(name,".sgrd$", sep=''),
                             full.names = T, recursive = T)

        print(sprintf("Using %s as input DEM", in.DEM))

      }else if (length(name) > 1){ # need to mosaic and transform
        print('Multiple grids Necessary...')
        original.DEM <- GetTilePathsHS(name, in.DEM)
        in.DEM <- MosaicAndWarp(gridnames = name, DEM.path = in.DEM,
                                saga.env = saga.env, outputCRS = projected.CRS)
      }

    }else if (toupper(DEM.source) %in% c('CDED','NED', 'CDEM', 'CDSM')){
      print(DEM.source) #debug only
      in.DEM <- OverlayDEM(point, DEM.dir=DEM.path, output.dir=saga.env$workspace,
                           product = DEM.source, tol=dem.clip.square) # clip from DEM

      dem.clip.square <- 0
    }else{
      print("Could not find suitable DEM in folder")
      return()
    }
  }

  # Test if point is outside of raster limits (e.g. in USA)
  value <- SampleRasterRS(point = point, grid = in.DEM, saga.env=saga.env )

  nodata <- F
  if (is.na(value)){
    print("station outside of DEM coverage, switching to NED data")
    in.DEM <- OverlayDEM(point, DEM.dir=DEM.path, output.dir=saga.env$workspace,
                         product = 'NED', tol=dem.clip.square)
    DEM.source <- "NED"
  }

  # Clip grid to point
  if (dem.clip.square > 0){
    ClipGridRS(point, in.DEM, 'clipped.sgrd', saga.env, tol=dem.clip.square, ...)
  }

  # Fill Sinks
  FillSinksRS(ifelse((dem.clip.square > 0),'clipped.sgrd',in.DEM),
              'filled.sgrd', saga.env, MINSLOPE=0.01, ...)

  # Upslope
  snapped.pt <- SnapToPourPoint(point, 'filled.sgrd', saga.env=saga.env,
                                upstream.threshold=upstream.threshold)

  point <- snapped.pt$station
  snap.dist <- snapped.pt$distance
  pointbuffer <- pb.init
  # Point to Grid
  Point2GridRS(point, 'filled.sgrd', 'target.sgrd',  saga.env,
               buffer.dist=pointbuffer, ...)

  upslope <- UpslopeAreaRS('filled.sgrd', 'target.sgrd', 'upslope.sgrd',
                           saga.env, method=method, ...)

  if (method==1){
    upslope <- GridThresholdRS(upslope, threshold.value = 0, saga.env=envi)
  }
  # Convert to polygon
  Grid2PolyRS(upslope, final.name, saga.env, ...)

  # calculate area
  area <- GridVolumeRS(upslope, level = 99, saga.env=envi)*1e-6

  # put workspace back the way it was
  setwd(oldwd)
  meta <- data.frame(
    station_number=station.name,
    final.name=final.name,
    area=area,
    pointbuffer=pointbuffer,
    nodata=nodata,
    snap.dist=snap.dist,
    DEM = DEM.source,
    cell_acc = upstream.threshold,
    method = upstr_method)
  write.csv(meta, gsub("shp", "csv", final.name), row.names = F, quote=F)
}

#===============================================================================
#' @title Create HydroBASINS catchment for a hydrometric station on a lake
#'
#' @description Create Hydrobasins upstream boundaries for a hydrometric station
#' on a lake using the HydroBASINS data product
#'
#' @param station_point a Spatialpointsdataframe of the hydrometric station
#'
#' @param saga.env an rsaga environment object
#'
#' @param output_folder where to save final output
#'
#' @return A list containing two spatialpolygons data frames (the lake and
#' the catchment area)
#'
#' @export
#===============================================================================
HYBASBasinLimits_Lake <- function(station_point, HYBAS, lakes_folder,
                                  output_folder, stn_num_col='station_number',
                                  stn_name_col='station_name'){

  station_name <- as.data.frame(station_point)[,stn_name_col]
  print(station_name)
  station_number <- as.data.frame(station_point)[,stn_num_col]
  waterbody <- ParseStationName(station_name)
  name <- paste0(waterbody,".*.shp")
  match <-list.files(lakes_folder, pattern=name, full.names=T)

  # get shapefile if it exists, otherwise download & convert kml
  if(length(match)==0){
    lac_geom <- CGNDBHydroKML(waterbody,
                              file.path(lakes_folder,gsub("\\*\\.shp$", "kml",name)))
  }else if (length(match)>=1){
    if (length(match)>1){
      match <-list.files(lakes_folder, pattern=paste0(waterbody,"_polygon.shp"),
                         full.names=T)[1]
      warning(sprintf("Multiple matching lake shapefiles found for %s.  Using %s",
                      station_name, match))
    }
    lac_geom <- readOGR(match, verbose=F, stringsAsFactors = FALSE)
    if(nrow(lac_geom)==0){return(NULL)}
    lac_geom <- outerRing(lac_geom)
}

  # reproject into planar coords
  albers <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0
                 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

  lac_geom_p <- sp::spTransform(lac_geom, sp::CRS(albers))
  station_point_p <- sp::spTransform(station_point, sp::CRS(albers))

  # find nearest lake polygon
  distance <- gDistance(station_point_p, lac_geom_p, byid=TRUE)
  min_dist <- min(distance)
  print(sprintf("Lake-station distance is %.01f m", min_dist))
  lac_geom <- lac_geom[which.min(distance),]

  # get intersect with hydrobasins
  basins <- sp::over(lac_geom, HYBAS, returnList = T)[[1]]

  # Find 'outlet' hybas and return all upstream
  outlet <- basins[which.max(basins$UP_AREA),]
  watershed <- FindAllUpstreamSubBasins(HYBAS, outlet$HYBAS_ID)

  #also find stn containing hybas (in case lake polygon misses it somehow)
  containing <- sp::over(station_point, HYBAS, returnList = T)[[1]]
  stn_up <- FindAllUpstreamSubBasins(HYBAS, containing$HYBAS_ID)

  hyb_select <- union(c(outlet$HYBAS_ID, watershed), c(stn_up, containing$HYBAS_ID))
  HYB_out <- HYBAS[HYBAS$HYBAS_ID %in% hyb_select,]
  print(outlet$HYBAS_ID)


  if (!missing(output_folder)){
    out_file <- file.path(output_folder, paste0(station_number, "_HYBAS_Basin.shp"))
    writeOGR(obj = HYB_out, dsn = out_file, layer = basename(out_file) ,
             driver = 'ESRI Shapefile')
    return(list(lake=lac_geom, hybas=HYB_out))
  }else{
    return(list(lake=lac_geom, hybas=HYB_out))
  }
}




#===============================================================================
#' @title Calculate basin area for a hydrometric station on a lake using a DEM
#'
#' @description Create a drainage basin shapefile using a DEM and a
#'
#' @param point a SpatialPointsDataFrame corresponding to a hydrometric station.
#'  Must be in same coordinates system as DEM and coordinate system
#' must be projected (not lat/long).  Must have longitude and latitude attributes.
#'
#' @param saga.env an rsaga environment object
#'
#' @param outdir Directory to output final upslope shapefile
#'
#' @param DEM.source character, one of: c('CDED', 'NED', 'CDEM, CDSM', 'SHEDS').
#'  Ignored if a DEM file is supplied to DEM.path
#'
#' @param DEM.path One of:
#'  (1) a file path to a directory containing DEM files either in the the format
#'   n%%w0%%_con_grid.sgrd (if DEM.source is 'SHEDS') or a directory to which DEM
#'   tiles will be downloaded (if DEM.source='CDEM', 'CDED' etc.),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a
#'  projected coordinate system and the coordinate system should match that of
#'  the point (e.g. Canada Albers Conformal Conic)
#'
#'
#' @return A l
#'
#' @export
#===============================================================================
LakeDEM <- function(station_point, lake_poly, lake_hybas_limit, DEM.path, outdir,
                    saga.env,
                    DEM.source='SHEDS',
                    stn_name_col = 'station_name',
                    stn_num_col = 'station_number'){
  # set params
  station_name   <- station_point@data[, stn_name_col]
  print(station_name)
  station_number <- station_point@data[, stn_num_col]
  albers <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0
                 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

  final.name <- file.path(outdir, paste(station_number,"_upslope.shp",sep=''))

  # get snap distance for output
  snap_dist <- as.character(
    gDistance(sp::spTransform(station_point, sp::CRS(albers)),
                         sp::spTransform(lake_poly, sp::CRS(albers)), byid=TRUE))

  # find most downstream hybas that overlaps lake

  lake_overlap <- over(lake_poly, lake_hybas_limit, returnlist=T)
  if (class(lake_overlap)=='list'){lake_overlap <- lake_overlap[[1]]}
  lowest_id <- lake_overlap$HYBAS_ID[which.max(lake_overlap$UP_AREA)]

  #hyb_max <- which.max(lake_hybas_limit$UP_AREA)
  outlet <- lake_hybas_limit[lake_hybas_limit$HYBAS_ID==lowest_id,]

  # clip out front of lake
  lake_clip <- gIntersection(outlet, lake_poly)
  lake_clip <- SpatialPolygonsDataFrame(lake_clip, data=data.frame(what='lake_clip'))
  lake_clip_p <- sp::spTransform(lake_clip, sp::CRS(albers)) #reproject

  # Make DEM
  point <- gCentroid(lake_clip)
  point <- SpatialPointsDataFrame(point@coords, data=data.frame(what='centroid'),
                                  proj4string = lake_clip@proj4string)
  point_p <- sp::spTransform(point, sp::CRS(albers)) #reproject
  outlet_p <- sp::spTransform(outlet, sp::CRS(albers))

  dem.clip.square <- max(outlet_p@bbox[,2] - outlet_p@bbox[,1])*1.25

  in.DEM <- gsub("\\.sdat$", "\\.sgrd", DEM.path)

  if (tolower(DEM.source)=='sheds'){
    name <- HydroMosaic(point@coords[,1], point@coords[,2], tol=dem.clip.square)
    if (length(name)==1){
      check.if.exists <- list.files(in.DEM, pattern=paste(name, ".sgrd", sep=''),
                                    full.names = T, recursive = T)
      if (length(check.if.exists) ==0){  #  if sgrid doesn't exist, create it
        original.DEM <- GetTilePathsHS(name, in.DEM)
        print(sprintf("Converting %s DEM to sgrd...", original.DEM))
        dstfile <- gsub("_con$","_con\\.sdat",original.DEM)
        gdal_warp2SAGA(original.DEM, outputCRS = projected.CRS, dstfile = dstfile)

      }
      in.DEM <- list.files(in.DEM, pattern=paste(name,".sgrd$", sep=''),
                           full.names = T, recursive = T)

      print(sprintf("Using %s as input DEM", in.DEM))

    }else if (length(name) > 1){ # need to mosaic and transform
      print('Multiple grids Necessary...')
      original.DEM <- GetTilePathsHS(name, in.DEM)
      in.DEM <- MosaicAndWarp(gridnames = name, DEM.path = in.DEM,
                              saga.env = saga.env, outputCRS = albers)
    }

  }else if (toupper(DEM.source) %in% c('CDED','NED', 'CDEM', 'CDSM')){
    print(DEM.source) #debug only
    in.DEM <- OverlayDEM(point, DEM.dir=DEM.path, output.dir=saga.env$workspace,
                         product = DEM.source, tol=dem.clip.square) # clip from DEM

  }else{
    print("Could not find suitable DEM in folder")
    return()
  }

  # Test if point is outside of raster limits (e.g. in USA)
  value <- SampleRasterRS(point = point_p, grid = in.DEM, saga.env=saga.env )
  nodata <- FALSE

  if (is.na(value)){
    print("station outside of CDED coverage, switching to NED data")
    nodata <- TRUE
    in.DEM <- OverlayDEM(point, DEM.dir=DEM.path, output.dir=saga.env$workspace,
                         product = 'NED', tol=dem.clip.square)
  }

  # Clip grid to point
  if (dem.clip.square > 0){
    ClipGridRS(point_p, in.DEM, 'clipped.sgrd', saga.env, tol=dem.clip.square)
  }

  # Fill Sinks
  FillSinksRS(ifelse((dem.clip.square > 0),'clipped.sgrd',in.DEM),
              'filled.sgrd', saga.env, MINSLOPE=0.01)

  # rasterize
  target <- Point2GridRS(lake_clip_p, 'filled.sgrd', 'target.sgrd',
                         saga.env, buffer.dist=0)

  upslope <- UpslopeAreaRS('filled.sgrd', 'target.sgrd', 'upslope.sgrd',
                           saga.env, method=0)

  # Convert to polygon (save final shape)
  Grid2PolyRS(upslope, final.name, saga.env)

  # calculate area
  area <- GridVolumeRS(upslope, level = 99, saga.env=envi)*1e-6

  meta <- data.frame(
    station_number=station_number,
    final.name = final.name,
    area = area,
    pointbuffer=0,
    nodata = nodata,
    snap.dist = snap_dist,
    DEM = DEM.source,
    cell_acc = 0,
    method = 'D8')

  # write output
  write.csv(meta, sub("shp$", "csv", final.name), row.names = F, quote=F)
}



#===============================================================================
#' @title Make upstream basin boundaries for a station on a lake
#'
#' @description This function acts a wrapper for two other functions
#'
#' @param point a SpatialPointsDataFrame corresponding to a hydrometric station.
#'  Must be in same coordinates system as DEM and coordinate system
#' must be projected (not lat/long).  Must have longitude and latitude attributes.
#'
#' @param saga.env an rsaga environment object
#'
#' @param outdir_HYBAS Directory to output final upslope hydrobasins boundaries
#'
#' @param DEM.source character, one of: c('CDED', 'NED', 'CDEM, CDSM', 'SHEDS').
#'  Ignored if a DEM file is supplied to DEM.path
#'
#' @param DEM.path One of:
#'  (1) a file path to a directory containing DEM files either in the the format
#'   n%%w0%%_con_grid.sgrd (if DEM.source is 'SHEDS') or a directory to which DEM
#'   tiles will be downloaded (if DEM.source='CDEM', 'CDED' etc.),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a
#'  projected coordinate system and the coordinate system should match that of
#'  the point (e.g. Canada Albers Conformal Conic)
#'
#'
#' @export
#===============================================================================
Lake_Station <-  function(station_point, HYBAS, lakes_folder, DEM.path,
                          DEM.source, outdir_DEM,
                          outdir_HYBAS, saga.env,
                          stn_num_col='station_number',
                          stn_name_col='station_name' ){

  shapes <-  HYBASBasinLimits_Lake (station_point=station_point,
                                    HYBAS=HYBAS,
                                    lakes_folder=lakes_folder,
                                    output_folder=outdir_HYBAS,
                                    stn_num_col=stn_num_col,
                                    stn_name_col=stn_name_col)

  if (is.null(shapes)){return(NULL)}

  LakeDEM(station = station_point,
          lake_poly = shapes$lake,
          lake_hybas_limit = shapes$hybas,
          DEM.path = DEM.path,
          DEM.source = DEM.source,
          outdir = outdir_DEM,
          saga.env = saga.env,
          stn_num_col=stn_num_col,
          stn_name_col=stn_name_col)

}



