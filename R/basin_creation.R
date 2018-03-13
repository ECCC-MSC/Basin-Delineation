
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
    next2up <- FindUpstreamSubBasins(HYBAS=HYBAS, HYBAS.ID = df[,hyb.id], ignore.endorheic = T)
    next2up <- HYBAS[HYBAS$HYBAS_ID %in% next2up,]
    nextup <- next2up@data[which.min(next2up$UP_AREA), "HYBAS_ID"] # find ID of next-up on smaller branch
    ids <- FindAllUpstreamSubBasins(HYBAS=HYBAS, HYBAS.ID = nextup)
    ids <- c(ids, nextup)

  }else if (code %in% c("Pcm")){
    next2up <- FindUpstreamSubBasins(HYBAS=HYBAS, HYBAS.ID = df[,hyb.id], ignore.endorheic = T)
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
  outname <- paste(df[,df.id], "HYBAS_Basin.shp", sep='_')
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


  # add in data to attributes table
  data <- DEM_data[,c('station_number', 'pointbuffer', 'nodata', 'snap.dist',
                      'DEM', 'cell_acc', 'method')]
  data$area <- round(raster::area(output)*1e-6,2)
  rownames(data) <- sapply(output@polygons, slot, name='ID')
  output <- SpatialPolygonsDataFrame(output, data)

  output <- sp::spTransform(output, CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

  # write output
  file_out <- file.path(out_folder, paste0(station,"_cntrb_basin.shp"))
  rgdal::writeOGR(obj = output,
           dsn = file_out,
           layer = basename(file_out),
           driver = 'ESRI Shapefile')
  write.csv(output@data, gsub("shp", "csv", file_out), row.names=F, quote=F)
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
      name <- HydroMosaic(point@data$longitude, point@data$latitude, tol=dem.clip.square)
      if (length(name)==1){
        check.if.exists <- list.files(in.DEM, pattern=paste(name, ".sgrd", sep=''), full.names = T, recursive = T)
        if (length(check.if.exists) ==0){  #  if sgrid doesn't exist, create it
          original.DEM <- GetTilePathsHS(name, in.DEM)
          print(sprintf("Converting %s DEM to sgrd...", original.DEM))
          dstfile <- gsub("_con$","_con\\.sdat",original.DEM)
          gdal_warp2SAGA(original.DEM, outputCRS = projected.CRS, dstfile = dstfile)

        }
        in.DEM <- list.files(in.DEM, pattern=paste(name,".sgrd$", sep=''), full.names = T, recursive = T)
        print(sprintf("Using %s as input DEM", in.DEM))

      }else if (length(name) > 1){ # need to mosaic and transform
        print('Multiple grids Necessary...')
        original.DEM <- GetTilePathsHS(name, in.DEM)
        in.DEM <- MosaicAndWarp(gridnames = name, DEM.path = in.DEM, saga.env = saga.env, outputCRS = projected.CRS)
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
    print("station outside of CDED coverage, switching to NED data")
    in.DEM <- OverlayDEM(point, DEM.dir=DEM.path, output.dir=saga.env$workspace,
                         product = 'NED', tol=dem.clip.square)
  }

  # Clip grid to point
  if (dem.clip.square > 0){
    ClipGridRS(point, in.DEM, 'clipped.sgrd', saga.env, tol=dem.clip.square, ...)
  }

  # Fill Sinks
  FillSinksRS(ifelse((dem.clip.square > 0),'clipped.sgrd',in.DEM),
              'filled.sgrd', saga.env, MINSLOPE=0.01, ...)

  # Upslope
  snapped.pt <- SnapToPourPoint(point, 'filled.sgrd', saga.env=saga.env, upstream.threshold=upstream.threshold)
  point <- snapped.pt$station
  snap.dist <- snapped.pt$distance
  pointbuffer <- pb.init
  # Point to Grid
  Point2GridRS(point, 'filled.sgrd', 'target.sgrd',  saga.env, buffer.dist=pointbuffer, ...)

  upslope <- UpslopeAreaRS('filled.sgrd', 'target.sgrd', 'upslope.sgrd', saga.env, method=method, ...)
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






