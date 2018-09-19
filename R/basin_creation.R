
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
#' @param hyb_id character, name of column containing ID of HydroBASINS sub-basin
#'  that contains the station
#'
#' @param stn_id character, name of column containing the station number
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
         hyb_id="HYBAS_ID", stn_id = "station_number"){
  print(paste(df[, stn_id], code))

  # get upstream HYBAS
  if (code %in% c("Pcb", "Ic")){
    ids <- FindAllUpstreamSubBasins(HYBAS = HYBAS, HYBAS.ID = df[, hyb_id])

  }else if (code %in% c("Pct")){
    next2up <- FindUpstreamSubBasins(HYBAS = HYBAS,
                                     HYBAS.ID = df[, hyb_id],
                                     ignore.endorheic = T)

    next2up <- HYBAS[HYBAS$HYBAS_ID %in% next2up, ]
    nextup <- next2up@data[which.min(next2up$UP_AREA), "HYBAS_ID"] # find ID of next-up on smaller branch
    ids <- FindAllUpstreamSubBasins(HYBAS = HYBAS, HYBAS.ID = nextup)
    ids <- c(ids, nextup)

  }else if (code %in% c("Pcm")){
    next2up <- FindUpstreamSubBasins(HYBAS = HYBAS, HYBAS.ID = df[, hyb_id],
                                     ignore.endorheic = T)

    next2up <- HYBAS[HYBAS$HYBAS_ID %in% next2up, ]
    nextup <- next2up@data[which.max(next2up$UP_AREA), "HYBAS_ID"] # find ID of next-up on main branch
    ids <- FindAllUpstreamSubBasins(HYBAS = HYBAS, HYBAS.ID = nextup)
    ids <- c(ids, nextup)

  }else if (grepl("[Ts]", code) & !grepl("du", code)){
    ids <- character()

  }else{
    print("not a valid code")
    return()
  }
  ids <- c(df[, hyb_id], ids)

  # save output file
  out <- HYBAS[HYBAS$HYBAS_ID %in% ids, ]
  outname <- paste(df[, stn_id], "HYBAS_Basin.shp", sep = "_")
  outfile <- file.path(outdir, outname)
  writeOGR(obj = out, dsn = outfile, layer = gsub("\\.shp", "", outname),
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
#' @param station_point spatialpointsdataframe of station
#'
#' @export
#===============================================================================
CombineWatersheds <- function(station, shedspoly_folder, dempoly_folder,
                              out_folder, station_point){
  if (!class(station) == "character"){stop("station should be a character string")}

  # find hydrobasins polygon in folder
  HYB_shp <- list.files(shedspoly_folder,
                        pattern = paste0(station, ".*\\.shp$"), full.names = T)
  if (length(HYB_shp) > 1){warning("Multiple hydrobasins polygons found.
                                    Using the first one.")}

  # load and project hydrobasins polygon
  HYB <- rgdal::readOGR(HYB_shp[1])
  HYB <- sp::spTransform(HYB, sp::CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=40
                                       +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80
                                       +datum=NAD83 +units=m +no_defs"))
  # find DEM polygon in folder
  DEM_shp <- list.files(dempoly_folder, pattern = paste0(station, ".*\\.shp$"),
                        full.names = T)

  if (length(DEM_shp) == 0){
    file_out <- file.path(out_folder, paste0(station, "_cntrb_basin_HYBonly.shp"))
    output <- rgeos::gUnaryUnion(HYB[-which.max(HYB$UP_AREA), ])
    output <- sp::SpatialPolygonsDataFrame(output,
                               data = data.frame(station_number = station,
                                               pointbuffer      = NA,
                                               nodata           = NA,
                                               snap_dist        = NA,
                                               DEM              = "None",
                                               cell_acc         = NA,
                                               method           = "HYB Only",
                                               area             = area(output)))

    warnings(sprintf("station %s has no DEM data.  HYBAS data used", station))

  }else{
    if (length(DEM_shp) > 1){warnings("Multiple DEM polygons found. Using the
                                  first one.")}

    DEM <- rgdal::readOGR(DEM_shp[1])
    DEM_data_csv <- list.files(dempoly_folder, pattern = paste0(station, ".*\\.csv$"),
                               full.names = T)
    DEM_data <- read.csv(DEM_data_csv)

    # reproject DEM polygon
    DEM <- sp::spTransform(DEM, CRSobj = sp::CRS(HYB@proj4string@projargs))

    # clean shapes
    DEM <- zeroBuffer(DEM)
    HYB <- zeroBuffer(HYB)

    # Perform set operations
    A <- HYB[-which.max(HYB$UP_AREA), ] # all but most downstream
    C <- HYB[which.max(HYB$UP_AREA), ]  # most downstream

    # remove gaps
    if (!missing(station_point)){
      DEM <- fill_upstream_gaps(station = station_point, basin = DEM,
                                hybas = C)
    }

    # Buffer 1mm to clean up any weird slices
    DEM <- rgeos::gBuffer(DEM, width = 1e-3)

    # clip DEM basin to hydrobasin boudary
    output <- rgeos::gIntersection(DEM, HYB)

    if (nrow(A) != 0){output <- rgeos::gUnion(output, A)}

    # remove holes (inner rings) prior to disaggregation
    output <- outerRing(output)

    # remove 'floating' polygons
    output <- sp::disaggregate(output)
    output <- output[which.max(lapply(output@polygons, slot, name = "area")), ]


    #  remove gaps in most downstream HYBAS
    # if (nrow(A) == 0 & !missing(station_point)){
    #   DEM <- fill_upstream_gaps(station = station_point, basin = output,
    #                                hybas = C)
    #
    # }else if (nrow(A) != 0 & !missing(station_point)){
    #   output <- fill_upstream_gaps(station = station_point, basin = output,
    #                                hybas = HYB, additive = F)
    #
    # }

    # add in data to attributes table
    snap <- ifelse('snap_dist' %in% names(DEM_data), 'snap_dist', 'snap.dist')
    data <- DEM_data[, c("station_number", "pointbuffer", "nodata", snap,
                        "DEM", "cell_acc", "method")]
    data$area <- round(raster::area(output) * 1e-6, 2)
    rownames(data) <- sapply(output@polygons, slot, name = "ID")
    output <- sp::SpatialPolygonsDataFrame(output, data)

    output <- sp::spTransform(output,
              CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

    file_out <- file.path(out_folder, paste0(station, "_cntrb_basin.shp"))
  }

  # write output
  rgdal::writeOGR(obj = output,
           dsn = file_out,
           layer = basename(file_out),
           driver = "ESRI Shapefile")

  write.csv(output@data, sub("shp$", "csv", file_out), row.names = F, quote = F)

  print("done")
  }

#===============================================================================
#' @title Calculate basin area using a DEM
#'
#' @description Create a drainage basin shapefile using a DEM and a
#' @param point a SpatialPointsDataFrame corresponding to a hydrometric station.
#'  Must be in same coordinates system as DEM and coordinate system
#' must be projected (not lat/long).  Must have longitude and latitude attributes.
#'
#' @param saga_env an rsaga.environment object
#'
#' @param outdir Directory to output final upslope shapefile
#'
#' @param pointbuffer numeric, how much should the point be buffered (calculates
#' upslope basin of buffered area)
#'
#' @param DEM_source character, one of: c('CDED', 'NED', 'CDEM, CDSM', 'SHEDS').
#'  Ignored if a DEM file is supplied to DEM_path
#'
#' @param DEM_path One of:
#'  (1) a file path to a directory containing DEM files either in the the format
#'   n%%w0%%_con_grid.sgrd (if DEM_source is 'SHEDS') or a directory to which DEM
#'   tiles will be downloaded (if DEM_source='CDEM', 'CDED' etc.),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a
#'  projected coordinate system and the coordinate system should match that of
#'  the point (e.g. Canada Albers Conformal Conic)
#'
#' @param clip_size How big (in m) a clip should be generated from the
#' original DEM (too big doesn't work and is slower)
#'
#' @param projected_CRS character vector (proj4) specifying "working" CRS in
#' which area and distance calculations are to be done
#'
#' @param method character, one of ('D8', 'DINF'), specifying the hydrological
#' pathing model
#'
#' @param upstr_threshold integer, number of upstream cells required to
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
DEMDrainageBasin <- function(point, DEM_path, DEM_source="NTS", saga_env,
                             outdir, pointbuffer=50, clip_size=50000,
                             projected_CRS,
                             upstr_method="D8", upstr_threshold=35000,
                             prelim_basin=NULL, ...){

  # set some parameters
  method <- switch(toupper(upstr_method), "D8" = 0, "DINF" = 1)
  snapped <- FALSE
  pb_init <- pointbuffer

  if (!is.null(prelim_basin)){
    clipping_layer <- prelim_basin
  }else{
    clipping_layer <- point
  }

  #sanity checks
  if (missing(projected_CRS)){
    projected_CRS <- GetProj4("AlbersEqualAreaConic")
  }
  in_DEM <- gsub("\\.sdat$", "\\.sgrd", DEM_path)

  # prepare workspace
  oldwd <- getwd()
  setwd(saga_env$workspace)
  outdir <- gsub("[\\/]$", "", outdir)

  # get name of station
  station.name <- point@data[1, grepl("^\\d{2}[[:alpha:]]{2}\\d{3}$", point@data)]
  print(station.name)
  final_name <- file.path(outdir, paste(station.name, "_upslope.shp", sep = ""))

  #if point not in proper crs, convert it
  if (point@proj4string@projargs != projected_CRS){
    point <- sp::spTransform(point, CRSobj = sp::CRS(projected_CRS))
  }

  # Find missing DEM if necessary / Convert to sgrd if necessary.
  if (!grepl("\\.sgrd$", in_DEM)){  # if not an SGRD file (then we expect a directory or a list of names)

    if (tolower(DEM_source) == "sheds"){
      name <- HydroMosaic(point@data$longitude, point@data$latitude,
                          tol = clip_size)

      if (length(name) == 1){
        matching_files <- list.files(in_DEM,
                                      pattern = paste(name, ".sgrd", sep = ""),
                                      full.names = T, recursive = T)

        if (length(matching_files) == 0){  #  if sgrid doesn't exist, create it
          original_DEM <- GetTilePathsHS(name, in_DEM)
          print(sprintf("Converting %s DEM to sgrd...", original_DEM))
          dstfile <- gsub("_con$", "_con\\.sdat", original_DEM)
          gdal_warp2SAGA(original_DEM, output_crs = projected_CRS, dstfile = dstfile)
        }

        in_DEM <- list.files(in_DEM, pattern = paste(name, ".sgrd$", sep = ""),
                             full.names = T, recursive = T)

        print(sprintf("Using %s as input DEM", in_DEM))

      }else if (length(name) > 1){ # need to mosaic and transform
        print("Multiple grids Necessary...")

        in_DEM <- MosaicAndWarp(gridnames = name, DEM_path = in_DEM,
                                saga_env = saga_env, output_crs = projected_CRS)
      }

    }else if (toupper(DEM_source) %in% c("CDED", "NED", "CDEM", "CDSM")){

      in_DEM <- OverlayDEM(clipping_layer, DEM_dir = DEM_path,
                           output_dir = saga_env$workspace, product = DEM_source,
                           tol = clip_size) # clip from DEM

      if (DEM_source != "NED"){clip_size <- 0}
    }else{
      print("Could not find suitable DEM in folder")
      return()
    }
  }

  # Clip grid to point
  if (clip_size > 0  ){
    ClipGridRS(clipping_layer, in_DEM, "clipped.sgrd", saga_env,
               tol = clip_size, ...)
  }

  # Fill Sinks
  FillSinksRS(in_DEM   = ifelse( (clip_size > 0), "clipped.sgrd", in_DEM),
              out_DEM  = "filled.sgrd",
              saga_env = saga_env, MINSLOPE = 0.01, ...)

  # Upslope
  snapped_pt <- SnapToPourPoint(point, "filled.sgrd", saga_env = saga_env,
                                upstr_threshold = upstr_threshold)

  point <- snapped_pt$station
  snap_dist <- snapped_pt$distance
  pointbuffer <- pb_init

  # Point to Grid
  Point2GridRS(point, "filled.sgrd", "target.sgrd",  saga_env,
               buffer_dist = pointbuffer, ...)

  upslope <- UpslopeAreaRS("filled.sgrd", "target.sgrd", "upslope.sgrd",
                           saga_env, method = method, ...)

  if (method == 1){
    upslope <- GridThresholdRS(upslope, threshold_value = 0, saga_env = saga_env)
  }
  # Convert to polygon
  Grid2PolyRS(upslope, final_name, saga_env, ...)

  # calculate area
  #area <- GridVolumeRS(upslope, level = 99, saga_env = saga_env) * 1e-6 # changed to polygon area for speed
  basin <- rgdal::readOGR(final_name)
  area <- sum(area(basin)) * 1e-6

  # put workspace back the way it was
  setwd(oldwd)
  meta <- data.frame(
    station_number = station.name,
    final_name     = final_name,
    area           = area,
    pointbuffer    = pointbuffer,
    nodata         = NA,
    snap_dist      = snap_dist,
    DEM            = DEM_source,
    cell_acc       = upstr_threshold,
    method         = upstr_method)

  write.csv(meta, gsub("shp", "csv", final_name), row.names = F, quote = F)
}

#===============================================================================
#' @title Create HydroBASINS catchment for a hydrometric station on a lake
#'
#' @description Create Hydrobasins upstream boundaries for a hydrometric station
#' on a lake using the HydroBASINS data product
#'
#' @param station_point a Spatialpointsdataframe of the hydrometric station
#'
#' @param saga_env an rsaga.environment object
#'
#' @param output_folder where to save final output
#'
#' @return A list containing two spatialpolygons data frames (the lake and
#' the catchment area)
#'
#' @export
#===============================================================================
HYBASBasinLimits_Lake <- function(station_point, HYBAS, lakes_folder,
                                  output_folder, code,
                                  stn_num_col="station_number",
                                  stn_name_col="station_name"){

  station_name <- as.data.frame(station_point)[, stn_name_col]
  print(station_name)

  station_number <- as.data.frame(station_point)[, stn_num_col]
  waterbody <- ParseStationName(station_name)
  name <- paste0(waterbody, ".*.shp")
  match <- list.files(lakes_folder, pattern = name, full.names = T)

  # get shapefile if it exists, otherwise download & convert kml
  if (length(match) == 0){
    lac_geom <- CGNDBHydroKML(waterbody,
                              file.path(lakes_folder,
                                        gsub("\\*\\.shp$", "kml", name)))

  }else if (length(match) >= 1){
    if (length(match) > 1){
      match <- list.files(lakes_folder, pattern = paste0(waterbody, "_polygon.shp"),
                         full.names = T)[1]
      warning(sprintf("Multiple matching lake shapefiles found for %s.  Using %s",
                      station_name, match))
    }

    lac_geom <- rgdal::readOGR(match, verbose=F, stringsAsFactors = FALSE)

    if (nrow(lac_geom) == 0){return(NULL)}

    lac_geom <- outerRing(lac_geom)
}

  # reproject into planar coords
  albers <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0
                 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

  lac_geom_p <- sp::spTransform(lac_geom, sp::CRS(albers))
  station_point_p <- sp::spTransform(station_point, sp::CRS(albers))

  # find nearest lake polygon
  distance <- rgeos::gDistance(station_point_p, lac_geom_p, byid = TRUE)
  min_dist <- min(distance)
  print(sprintf("Lake-station distance is %.01f m", min_dist))
  lac_geom <- lac_geom[which.min(distance), ]

  # get intersect with hydrobasins
  basins <- sp::over(lac_geom, HYBAS, returnList = T)[[1]]

  # for 's' coded stations, return containing hybas ONLY
  if (grepl('s', code)){
    HYB_out <- HYBAS[HYBAS$HYBAS_ID == basins$HYBAS_ID[1], ]
    return(list(lake = lac_geom, hybas = HYB_out))
  }

  # Find 'outlet' hybas and return all upstream
  outlet <- basins[which.max(basins$UP_AREA), ]
  watershed <- FindAllUpstreamSubBasins(HYBAS, outlet$HYBAS_ID)

  #also find stn containing hybas (in case lake polygon misses it somehow)
  containing <- sp::over(station_point, HYBAS, returnList = T)[[1]]
  stn_up <- FindAllUpstreamSubBasins(HYBAS, containing$HYBAS_ID)

  hyb_select <- union(c(outlet$HYBAS_ID, watershed), c(stn_up, containing$HYBAS_ID))
  HYB_out <- HYBAS[HYBAS$HYBAS_ID %in% hyb_select, ]
  print(outlet$HYBAS_ID)


  if (!missing(output_folder)){
    out_file <- file.path(output_folder, paste0(station_number, "_HYBAS_Basin.shp"))
    rgdal::writeOGR(obj = HYB_out, dsn = out_file, layer = basename(out_file),
             driver = "ESRI Shapefile")

    return(list(lake = lac_geom, hybas = HYB_out))

  }else{

    return(list(lake = lac_geom, hybas = HYB_out))

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
#' @param saga_env an rsaga.environment object
#'
#' @param outdir Directory to output final upslope shapefile
#'
#' @param DEM_source character, one of: c('CDED', 'NED', 'CDEM, CDSM', 'SHEDS').
#'  Ignored if a DEM file is supplied to DEM_path
#'
#' @param DEM_path One of:
#'  (1) a file path to a directory containing DEM files either in the the format
#'   n%%w0%%_con_grid.sgrd (if DEM_source is 'SHEDS') or a directory to which DEM
#'   tiles will be downloaded (if DEM_source='CDEM', 'CDED' etc.),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a
#'  projected coordinate system and the coordinate system should match that of
#'  the point (e.g. Canada Albers Conformal Conic)
#'
#'
#' @return A l
#'
#' @export
#===============================================================================
LakeDEM <- function(station_point, lake_poly, lake_hybas_limit, DEM_path, outdir,
                    saga_env,
                    DEM_source="SHEDS",
                    stn_name_col = "station_name",
                    stn_num_col = "station_number"){
  # set params
  station_name   <- station_point@data[, stn_name_col]
  print(station_name)
  station_number <- station_point@data[, stn_num_col]
  albers <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0
                 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

  final_name <- file.path(outdir, paste(station_number, "_upslope.shp", sep = ""))

  # get snap distance for output
  snap_dist <- as.character(
    gDistance(sp::spTransform(station_point, sp::CRS(albers)),
                         sp::spTransform(lake_poly, sp::CRS(albers)), byid = TRUE))

  # find most downstream hybas that overlaps lake

  lake_overlap <- sp::over(lake_poly, lake_hybas_limit, returnlist = T)
  if (class(lake_overlap) == "list"){lake_overlap <- lake_overlap[[1]]}
  lowest_id <- lake_overlap$HYBAS_ID[which.max(lake_overlap$UP_AREA)]

  #hyb_max <- which.max(lake_hybas_limit$UP_AREA)
  outlet <- lake_hybas_limit[lake_hybas_limit$HYBAS_ID == lowest_id, ]

  # clip out front of lake
  lake_clip <- rgeos::gIntersection(outlet, lake_poly)
  lake_clip <- sp::SpatialPolygonsDataFrame(lake_clip,
                                            data = data.frame(what = "lake_clip"))

  lake_clip_p <- sp::spTransform(lake_clip, sp::CRS(albers)) #reproject

  # Make DEM
  point <- rgeos::gCentroid(lake_clip)
  point <- sp::SpatialPointsDataFrame(point@coords, data = data.frame(what = "centroid"),
                                  proj4string = lake_clip@proj4string)

  point_p <- sp::spTransform(point, sp::CRS(albers)) #reproject
  outlet_p <- sp::spTransform(outlet, sp::CRS(albers))

  clip_size <- max(outlet_p@bbox[, 2] - outlet_p@bbox[, 1]) * 1.25

  in_DEM <- gsub("\\.sdat$", "\\.sgrd", DEM_path)

  if (tolower(DEM_source) == "sheds"){
    name <- HydroMosaic(point@coords[, 1], point@coords[, 2], tol = clip_size)

    if (length(name) == 1){
      matching_files <- list.files(in_DEM, pattern = paste(name, ".sgrd", sep = ""),
                                    full.names = T, recursive = T)
      if (length(matching_files) == 0){  #  if sgrid doesn't exist, create it
        original_DEM <- GetTilePathsHS(name, in_DEM)
        print(sprintf("Converting %s DEM to sgrd...", original_DEM))
        dstfile <- gsub("_con$", "_con\\.sdat", original_DEM)
        gdal_warp2SAGA(original_DEM, output_crs = projected_CRS, dstfile = dstfile)
      }

      in_DEM <- list.files(in_DEM, pattern = paste(name, ".sgrd$", sep = ""),
                           full.names = T, recursive = T)

      print(sprintf("Using %s as input DEM", in_DEM))

    }else if (length(name) > 1){ # need to mosaic and transform
      print("Multiple grids Necessary...")
      original_DEM <- GetTilePathsHS(name, in_DEM)
      in_DEM <- MosaicAndWarp(gridnames = name, DEM_path = in_DEM,
                              saga_env = saga_env, output_crs = albers)
    }

  }else if (toupper(DEM_source) %in% c("CDED", "NED", "CDEM", "CDSM")){
    print(DEM_source) #debug only
    in_DEM <- OverlayDEM(point, DEM_dir = DEM_path,
                         output_dir = saga_env$workspace,
                         product = DEM_source,
                         tol = clip_size) # clip from DEM

  }else{
    print("Could not find suitable DEM in folder")
    return()
  }

  # Test if point is outside of raster limits (e.g. in USA)
  value <- SampleRasterRS(point = point_p, grid = in_DEM, saga_env = saga_env )
  nodata <- FALSE

  if (is.na(value)){
    print("station outside of CDED coverage, switching to NED data")
    nodata <- TRUE
    in_DEM <- OverlayDEM(point, DEM_dir = DEM_path,
                         output_dir = saga_env$workspace,
                         product = "NED", tol = clip_size)
  }

  # Clip grid to point
  if (clip_size > 0){
    ClipGridRS(point_p, in_DEM, "clipped.sgrd", saga_env, tol = clip_size)
  }

  # Fill Sinks
  FillSinksRS(ifelse( (clip_size > 0), "clipped.sgrd", in_DEM),
              "filled.sgrd", saga_env, MINSLOPE = 0.01)

  # rasterize
  target <- Point2GridRS(lake_clip_p, "filled.sgrd", "target.sgrd",
                         saga_env, buffer_dist = 0)

  upslope <- UpslopeAreaRS("filled.sgrd", "target.sgrd", "upslope.sgrd",
                           saga_env, method = 0)

  # Convert to polygon (save final shape)
  Grid2PolyRS(upslope, final_name, saga_env)

  # calculate area
  area <- GridVolumeRS(upslope, level = 99, saga_env = saga_env) * 1e-6

  meta <- data.frame(
    station_number = station_number,
    final_name     = final_name,
    area           = area,
    pointbuffer    = 0,
    nodata         = nodata,
    snap_dist      = snap_dist,
    DEM            = DEM_source,
    cell_acc       = 0,
    method         = "D8")

  # write output
  write.csv(meta, sub("shp$", "csv", final_name), row.names = F, quote = F)
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
#' @param saga_env an rsaga.environment object
#'
#' @param outdir_HYBAS Directory to output final upslope hydrobasins boundaries
#'
#' @param DEM_source character, one of: c('CDED', 'NED', 'CDEM, CDSM', 'SHEDS').
#'  Ignored if a DEM file is supplied to DEM_path
#'
#' @param DEM_path One of:
#'  (1) a file path to a directory containing DEM files either in the the format
#'   n%%w0%%_con_grid.sgrd (if DEM_source is 'SHEDS') or a directory to which DEM
#'   tiles will be downloaded (if DEM_source='CDEM', 'CDED' etc.),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a
#'  projected coordinate system and the coordinate system should match that of
#'  the point (e.g. Canada Albers Conformal Conic)
#'
#'
#' @export
#===============================================================================
Lake_Station <-  function(station_point, HYBAS, lakes_folder, DEM_path,
                          DEM_source, outdir_DEM,
                          outdir_HYBAS, saga_env, code,
                          stn_num_col="station_number",
                          stn_name_col="station_name" ){

  shapes <-  HYBASBasinLimits_Lake(station_point  = station_point,
                                    HYBAS         = HYBAS,
                                    lakes_folder  = lakes_folder,
                                    output_folder = outdir_HYBAS,
                                    stn_num_col   = stn_num_col,
                                    stn_name_col  = stn_name_col,
                                    code          = code)

  if (is.null(shapes)){return(NULL)}

  LakeDEM(station          = station_point,
          lake_poly        = shapes$lake,
          lake_hybas_limit = shapes$hybas,
          DEM_path         = DEM_path,
          DEM_source       = DEM_source,
          outdir           = outdir_DEM,
          saga_env         = saga_env,
          stn_num_col      = stn_num_col,
          stn_name_col     = stn_name_col)

}



#===============================================================================
#' Create basin delineation from DEM and HydroBASINS (deprecated)
#'
#' @description Gives information about a station (upstream area, upstream
#' stations)
#'
#' @param station Either a SpatialPointsDataFrame or the path to a
#' Hydat-exported table with desired points or a character vector of station
#'  names (in this case a connection must be specified)
#'
#' @param con (optional) A open database connection.  Used if station is a
#' character vector of station names
#'
#' @param DEM_source either
#'
#' @param DEM.index index shapefile (optional)
#'
#' @param outdir (optional) A directory to which the output shapefile and report
#'  will be saved.  If not specified,
#' the function will return an R object without saving data to disk.
#'
#' @param pourpoints spatial points data frame with
#'
#' @param NCA (optional) Either a character string path or a
#' SpatialPolygonsDataFrame.  This is a polygon of the non-contributing areas.
#'  If provided, calculates effective drainage area
#'
#' @return A list of SpatialPolygonDataFrames
#'
#' @export
#===============================================================================
DelineateBasin <- function(station, con, outdir, HYBAS, DEM_path, DEM.index,
                           DEM_source, saga_env, NCA, interactive=F, pourpoints,
                           default.to.HYB=F,
                           projected_CRS="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40
                           +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83
                           +units=m +no_defs", ...){

  station <- ReadStationInformation(station, con = con)

  if (!missing(pourpoints)){
    station <- getPourPoint(station, pourpoints)
  }
  #=============
  # Calculations
  #=============

  # Get upslope area from HYBAS
  HYB.Poly.individual <- UpstreamHYBAS(station, HYBAS)
  HYB.Poly <- rgeos::gUnaryUnion(HYB.Poly.individual,
                                 id = HYB.Poly.individual$ID) # combined area

  no.upstream.hyb <- (length(HYB.Poly.individual) == 1)

  # Get first two upstream tiles (should be the first one for each ID)
  first.up <- unlist(by(
    data = HYB.Poly.individual@data,
    INDICES = HYB.Poly.individual@data$ID,
    FUN = function(x)
      as.character(x[which.max(as.character(x$UP_AREA)), c("HYBAS_ID")])))

  first.up.HYB <- HYB.Poly.individual[HYB.Poly.individual$HYBAS_ID %in%
                                    first.up & HYB.Poly.individual$ENDO != 2, ]

  first.up.HYB.endo <- HYB.Poly.individual[HYB.Poly.individual$HYBAS_ID %in%
                                             first.up, ]

  # find any endorheic basins
  endo.HYB <- HYB.Poly.individual@data[HYB.Poly.individual$HYBAS_ID %in%
                                 first.up & HYB.Poly.individual$ENDO == 2, "ID"]

  ID.list <- sapply(slot(HYB.Poly, "polygons"), function(x) slot(x, "ID"))
  endo <- ID.list %in% endo.HYB

  disjoint <- rep(F, length(endo))
  not.unlikely <- rep(F, length(endo))

  # find necessary coverage for DEM TODO[NB]: Clean part this up: its a mess!
  if (!missing(DEM.index)){
    DEM_path <- OverlayDEM(rgeos::gUnaryUnion(first.up.HYB),
                           tileindex = DEM.index, saga_env$workspace,
                           saga_env$workspace, product = DEM_source,
                           tilename = "NTS_SNRC")
  }else{
    DEM_path <- OverlayDEM(geom1 = rgeos::gUnaryUnion(first.up.HYB),
                           DEM.dir = saga_env$workspace,
                           output_dir = saga_env$workspace,
                           product = DEM_source, tol = 15e3
    )
  }

  # Get upslope area from DEM
  DEMresult <- UpslopeDEM(station, DEM_path = DEM_path, DEM_source = DEM_source,
                          saga_env = saga_env,
                          outdir = saga_env$workspace,
                          projected_CRS = projected_CRS, ...)

  if (DEMresult$nodata == T){ # if the point has no data
    DEM.poly.p <- sp::spTransform(HYB.Poly, CRS = sp::CRS(projected_CRS)) #use different CRS?
  }
  DEM.Poly <- DEMresult$final_name
  pointbuffer <- DEMresult$pointbuffer
  iterate_to <- DEMresult$iterate_to
  iterate_thres <- DEMresult$iterate_thres

  DEM.Poly.p <- invisible(rgdal::readOGR(DEM.Poly))
  DEM.Poly.p <- DEM.Poly.p[which.max(raster::area(DEM.Poly.p)), ] # take the big one.  others will likely have bad geometries

  # Clip 'nose' of basin
  HP.p <- sp::spTransform(HYB.Poly,
                          CRSobj = sp::CRS(DEM.Poly.p@proj4string@projargs))

  local.drainage.p <- rgeos::gIntersection(HP.p[1, ],
                                           DEM.Poly.p,
                                           byid = F,
                                           drop_lower_td = TRUE) # DEM area within containing hybas


  #=====================
  #  Logical Steps
  #=====================

  if (!no.upstream.hyb){
    # if DEM is disjoint with upslope area, do not include it.
    disjoint <- rgeos::gDisjoint(HP.p, local.drainage.p, byid = T)

    # if calculated total upslope area is less than half of an upslope HYBAS region, remove that HYBAS region
    not.unlikely <- raster::area(DEM.Poly.p) >
      0.5 * raster::area(sp::spTransform(first.up.HYB.endo,
                                       CRS(DEM.Poly.p@proj4string@projargs)))
  }

  keep <- ID.list[(!disjoint | endo) & (not.unlikely | endo)]
  HP.out <- HP.p[keep[-which(keep == "A")], ] # maybe can delete this?

  # merge desired parts of HYBAS polygons with local drainage
  if (pointbuffer >= iterate_to & default.to.HYB){
    # in this case, just use the HYBAS
    output <- rgeos::gUnaryUnion(HP.p)
  }else if ( (length(HP.out) > 0 & !no.upstream.hyb) |
            (pointbuffer >= iterate_to & !default.to.HYB)){
    # HP.p.merged <- rgeos::gUnaryUnion(HP.out)
    # output <- rgeos::gUnion(local.drainage.p, HP.out)
    ContainingHydrobasins <- HP.p[keep] #include containing ("A") hydrobasin
    DEMwithinHB <- rgeos::gIntersection(DEM.Poly.p, ContainingHydrobasins)
    output <- rgeos::gUnion(HP.out, DEMwithinHB)
  }else{
    output <- rgeos::gUnaryUnion(local.drainage.p)
  }

  dr_ar <- raster::area(output) * 1e-6

  # clean geometry with a zero-width buffer and remove holes / small rings
  output <- rgeos::gBuffer(output, byid = TRUE, width = 0)
  output <- outerRing(output)

  # Calculate non-contributing areas
  if (!missing(NCA)){

  }else{
    drainage.effective <- NA
  }

  # create attribute table
  station@data <- AddMissingColumns(station@data,
                                    c("drainage_area_gross", "drainage_area_effect"))

  data <- as.data.frame(list(  # make output table
    stn_number = station@data$station_number,
    drainar_gr = dr_ar,
    drainar_ef = drainage.effective,  # add notes?
    stn_dr_gr  = station@data$drainage_area_gross,
    stn_dr_ef  = station@data$drainage_area_effect,
    trgtRdius  = pointbuffer,
    trgtRdMax  = iterate_to,
    iterThres  = iterate_thres,
    snapped    = DEMresult$snapped,
    snapDist   = DEMresult$snap_dist,
    DEM_source = DEM_source
  ))

  output <-  sp::spTransform(SpatialPolygonsDataFrame(output, data), # output shape file with same CRS as station
                             CRSobj = sp::CRS(station@proj4string@projargs))

  # write output shapefile (rgdal::writeOGR)
  name <- paste(station$station_number, "_basin.shp")
  rgdal::writeOGR(output, dsn = file.path(outdir, name), layer = name,
                  driver = "ESRI Shapefile")

  # write output table (optional?)
  tablename <- paste(station$station_number, "_basin.csv")
  write.csv(output@data, file.path(outdir, tablename), row.names = F)

  #clean up (optional?)
  # file.remove(DEM.Poly)

  return(list(
    includedHYB = HP.out,
    dem_local   = local.drainage.p,
    out         = output,
    orig.DEM    = DEM.Poly.p,
    allHYB      = HP.p))
  }
