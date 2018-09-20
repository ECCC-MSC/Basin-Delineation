# troubleshooting RSAGA:
# - run command in commandline
# - sometimes the variables names are different between versions of SAGA, (e.g. INTERPOL instead of RESAMPLING)

#===============================================================================
#' @title Snap to pour point
#'
#' @description Finds the nearest location on a stream segment to a hydrological
#' station and returns the pourpoint as a spatialpointsdataframe with the same
#'  attributes as the input station.
#'
#' @param station spatialpointsdataframe representing a hydrological station
#'
#' @param in_DEM character path pointing to hydrologically conditioned DEM.
#' Must be in projected
#' coordinate system and, at a minimum, have filled sinks.
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @param upstr_threshold integer, number of upstream cells required to
#' initiate channel growth. A smaller number will yield more channels and will
#'  snap to smaller streams. Larger numbers will produce only major channel
#'   branches. For a 50k DEM, numbers between 20,000 and 50,000 are appropriate
#'   to reproduce the main channel network without many small branches.
#'
#' @return a SpatialPointsDataFrame with the same data as the input station,
#' but located on a DEM stream branch.
#'
#' @export
#===============================================================================
SnapToPourPoint <- function(station, in_DEM, saga_env, upstr_threshold = 35000){
  chann <- ChannelNetworkRS(in_DEM = in_DEM, upstr_threshold = upstr_threshold,
                            saga_env = envi, verbose = F)
  chann <- rgdal::readOGR(chann)
  newpoint <- SnapToNearest(station, chann)
  dist <- rgeos::gDistance(station, chann)

  return(list(station = newpoint, distance = dist))
}

#===============================================================================
#' @title Shapes to Grid
#'
#'
#' @description Turns a vector object into a grid ('grid_gridding' module 0)
#'
#' @param point A SpatialPoint or SpatialPointsDataFrame corresponding to a
#'  station
#'
#' @param grid.sys A SAGA *.sgrd file with the same grid system as the desired
#' output
#'
#' @param out_name Character string of output file name
#'
#' @param env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @param buffer_dist If equal to zero does not buffer.
#'
#' @keywords internal
#'
#' @export
#===============================================================================
Point2GridRS <- function(point, grid.sys, out_name,  saga_env, buffer_dist = 150,
                         verbose = F){
  print("Converting point to grid ...")
  # Write point as shapefile
  rgdal::writeOGR(point, dsn = file.path(saga_env$workspace, "target.shp"),
                  layer = "target.shp", driver = "ESRI Shapefile")
  # Convert shapefile to grid
  RSAGA::rsaga.geoprocessor(lib = "grid_gridding", module = 0, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    INPUT             = file.path(saga_env$workspace, "target.shp"),
    OUTPUT            = 0,
    GRID_TYPE         = 0,
    TARGET_DEFINITION = 1,  # Better to define by the other grid (less prone to error? but maybe slower??)
    TARGET_TEMPLATE   = grid.sys,
    GRID              = out_name
  ))

  if (buffer_dist > 0){
    # Grid buffer module 8
    RSAGA::rsaga.geoprocessor(lib = "grid_tools", module = 8, env = saga_env,
                              show.output.on.console = verbose,
                              check.module.exists = FALSE, param = list(
      FEATURES   = out_name,
      BUFFER     = out_name,
      DIST       = buffer_dist,
      BUFFERTYPE = 0
    ))
  }
}

#===============================================================================
#' @title Clip grid to object plus buffer (Clip Grids = 31)
#'
#' @param point Spatial* object
#'
#' @param in_DEM character path to tiled grid to clip
#'
#' @param out_DEM (optional) path to save output DEM, otherwise appends
#' '_clipped' to filename of in_DEM
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @param tol size of clipping box (metres) around point
#'
#' @keywords internal
#'
#' @export
#===============================================================================
ClipGridRS <- function(point, in_DEM, out_DEM, saga_env, tol, verbose=F){
  print("Clipping grid to vector ...")

  if (missing(out_DEM)) out_DEM <- gsub("\\.sgrd", "_clipped\\.sgrd", in_DEM)

  xmin <- point@bbox[1, 1]
  xmax <- point@bbox[1, 2]
  ymin <- point@bbox[2, 1]
  ymax <- point@bbox[2, 2]

  RSAGA::rsaga.geoprocessor(lib = "grid_tools", module = 31, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    GRIDS   = in_DEM,
    CLIPPED = out_DEM,
    XMIN    = xmin - tol,
    XMAX    = xmax + tol,
    YMIN    = ymin - tol,
    YMAX    = ymax + tol
  ))
}

#===============================================================================
#' @title Clip Polygon   (Clip Grids = 11)
#'
#' @param input SpatialPointsDataFrame
#'
#' @param clipping.layer character path to tiled grid to clip
#'
#' @param output (optional) path to save output DEM, otherwise appends '_clipped'
#' to filename of in_DEM
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @keywords internal
#'
#' @export
#===============================================================================
ClipPolygonRS <- function(input, clipping.layer, output, saga_env, verbose=F){
  print("Clipping Polygons...")
  if (missing(output)) output <- gsub("\\.sgrd", "_clipped\\.sgrd", input)

  RSAGA::rsaga.geoprocessor(lib = "shapes_polygons", module = 11, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
                              CLIP     = clipping.layer,
                              S_INPUT  = input,
                              S_OUTPUT = output,
                              DISSOLVE = 1
                              ))
}

#===============================================================================
#' @title Fill Sinks
#'
#'
#' @description Fill sinks according to Wang and Liu ('ta_preprocessor' module
#' 4 or 5 for XXL)
#'
#' @param in_DEM character path to grid to fill
#'
#' @param out_DEM (optional) path to save output DEM, otherwise appends
#' '_filled' to filename of in_DEM
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @param MINSLOPE minimum slope to maintain between cells
#'
#' @param XXL logical, whether or not to use the 'Fill sinks XXL' tool, suitable
#'  for large grids
#'
#' @export
#===============================================================================
FillSinksRS <- function(in_DEM, out_DEM, saga_env, MINSLOPE=0.01, verbose=F,
                        XXL=F){
  print("Filling sinks...")
  if (missing(out_DEM)) out_DEM <- gsub("\\.s[gd][ra][dt]", "_filled\\.sgrd",
                                        in_DEM)

  # Fill sinks (wang & Liu) = 4
  module <- ifelse(XXL, 5, 4)
  RSAGA::rsaga.geoprocessor(lib = "ta_preprocessor", module = module,
                            env = saga_env, show.output.on.console = verbose,
                    check.module.exists = FALSE, param = list(
                      ELEV     = in_DEM,
                      FILLED   = out_DEM,
                      MINSLOPE = MINSLOPE  # unsure how to choose optimal value
                            ))

  return(out_DEM)
}

#===============================================================================
#' @title Sink Removal
#'
#' @description Remove sinks by filling sinks or deepening drainage routes using
#' the method of Conrad (2001)
#'
#' @param in_DEM character path to grid to fill
#'
#' @param out_DEM (optional) path to save output DEM, otherwise appends '_filled'
#' to filename of in_DEM
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @param METHOD Available Choices: [0] Deepen Drainage Routes [1] Fill Sinks.
#'
#' @param THRESHOLD logical, whether or not to exclude deeper sinks from filling
#'
#' @param THRESHOLD.HEIGHT numeric, The parameter describes the maximum depth of
#' a sink to be considered for removal [map units]. This makes it possible to
#' exclude deeper sinks from filling.
#'
#' @keywords internal
#'
#' @export
#===============================================================================
BreachDepressionsRS <- function(in_DEM, out_DEM, saga_env, METHOD=0, THRESHOLD=0,
                                THRESHOLD.HEIGHT=100, verbose=F){
  print(paste(ifelse(METHOD, "Filling", "Breaching"), "sinks..."))

  if (missing(out_DEM)) out_DEM <- gsub("\\.sgrd", "_sinksremoved\\.sgrd", in_DEM)


  RSAGA::rsaga.geoprocessor(lib = "ta_preprocessor", module = 2, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
                              DEM         = in_DEM,
                              DEM_PREPROC = out_DEM,
                              METHOD      = METHOD,
                              THRESHOLD   = THRESHOLD,
                              THRSHEIGHT  = THRESHOLD.HEIGHT
                            ))
  return(out_DEM)
}

#===============================================================================
#' @title Upslope area
#'
#' @description Calculates upslope area and produces a grid ('ta_hydrology'
#' module 4)
#'
#' @param in_DEM Hydrologically appropriate DEM file (saga grid)
#'
#' @param target.grid saga grid
#'
#' @param out.GRD character string path to output grid
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @export
#===============================================================================
UpslopeAreaRS <- function(in_DEM, target.grid, out.GRD, saga_env, verbose=F,
                          method=1){
  print("Calculating upstream area...")
  RSAGA::rsaga.geoprocessor(lib = "ta_hydrology", module = 4, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    TARGET    = target.grid,
    ELEVATION = in_DEM,
    AREA      = out.GRD,
    METHOD    = method
  ))
  return(out.GRD)
}

#===============================================================================
#' @title Grid to polygon
#'
#' @description Creates a polygon from a grid ('shapes_grid' module 6:
#'  Vectorising Grid Classes)
#' @param in_grd character string path to SAGA grid file.  The grid should
#' @param out_poly character string path to output polygon shapefile
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param CLASS_ID numeric raster value that designates which grid cells to
#'  convert to polygon (cells with other values will be ignored)
#' @export
#===============================================================================
Grid2PolyRS <- function(in_grd, out_poly, saga_env, CLASS_ID=100, verbose=F){
  print("Converting grid to polygons...")
  RSAGA::rsaga.geoprocessor(lib = "shapes_grid", module = 6, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    GRID        = in_grd,
    POLYGONS    = out_poly,
    CLASS_ALL   = 0,
    CLASS_ID    = CLASS_ID,
    SPLIT       = 1,
    ALLVERTICES = 0
  ))
}

#===============================================================================
#' @title Difference grids
#' @description shapes_polygons' module 15
#===============================================================================
DifferenceRS <- function(shape, eraser, output, saga_env, verbose=F){
  print("Differencing shapes...")
  RSAGA::rsaga.geoprocessor(lib = "shapes_polygons", module = 15, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    A      = shape,
    B      = eraser,
    RESULT = output,
    SPLIT  = 0
    ))
}

#===============================================================================
#' @title Calculate Grid Volume
#' @description grid_calculus' module 2
#===============================================================================
GridVolumeRS <- function(grid, level, method=0, saga_env){
  print("Calculating Grid Volume...")

  volume <- invisible(capture.output(RSAGA::rsaga.geoprocessor(
    lib = "grid_calculus", module = 2, env = saga_env,
               check.module.exists = FALSE, param = list(
                              GRID   = grid,
                              METHOD = method,
                              LEVEL  = level
                            ))))

  volume <- as.numeric(gsub("Volume: ", "",
                            regmatches(volume, regexpr("^Volume: (\\d*\\.\\d*)",
                                                      volume))))
  return(volume)
}

#===============================================================================
#' @title Grid threshold
#===============================================================================
GridThresholdRS <- function(grid, threshold_value, method=0, saga_env, dstfile){
  print("Thresholding Grid...")

  if (missing(dstfile)) dstfile <- gsub("\\.(.{2,5})$", "_thres\\.sdat", grid)

  volume <- capture.output(RSAGA::rsaga.geoprocessor(lib = "grid_calculus",
                                                     module = 1, env = saga_env,
            check.module.exists = FALSE, param = list(
                   GRIDS   = grid,
                   FORMULA = sprintf("(g1>%d)*100", threshold_value),
                   TYPE    = 4,
                   RESULT  = dstfile
                    )))

  return(dstfile)
}
#===============================================================================
#' @title Mosaic
#' @description Builds mosaicked DEM using SAGA 'grid_tools' module 3
#' @param grids character vector of input grids
#' @param out_grid file path to output grid
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param xmin xmin
#' @param xmax xmax
#' @param ymin ymin
#' @param ymax ymax
#' @param cellsize cell size

#===============================================================================
MoisaicRS <- function(grids, out_grid, saga_env, xmin, xmax, ymin, ymax,
                      cellsize, verbose = F){
  print("Building Mosaic...")

  RSAGA::rsaga.geoprocessor(lib = "grid_tools", module = 3, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    GRIDS             = grids,
    TYPE              = 4,  # 2-byte signed integer
    INTERPOL          = 4,
    OVERLAP           = 0,
    MATCH             = 0,
    TARGET_DEFINITION = 0,
    TARGET_USER_XMIN  = xmin,
    TARGET_USER_XMAX  = xmax,
    TARGET_USER_YMIN  = ymin,
    TARGET_USER_YMAX  = ymax,
    TARGET_USER_SIZE  = cellsize,
    TARGET_OUT_GRID   = out_grid
    ))
}


#===============================================================================
#' @title Convert GRID to SAGA SGRD
#'
#' @description Convertes a GRID DEM into SAGA .sgrd/.sdat format.  Also
#' converts the CRS into a projected system.
#'
#' @param srcfile path to GRID folder (no trailing path separators)
#'
#' @param dstfile (optional) path to output SAGA grid file.  If not specified,
#' saves to same directory as input.  Note that the file name must have the
#' "sdat" file extension (not sgrid)
#'
#' @param s_srs source coordinate system.  Defaults to EPSG:4326 (WGS 84)
#'
#' @param output_crs character string specifying the coordinate system of the
#' output grid.  Examples include: "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40
#' +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#'
#' @keywords internal
#' @export
#===============================================================================
gdal_warp2SAGA <- function(srcfile, dstfile, output_crs, srcnodata = -32768,
                           dstnodata = -32767, s_srs = "EPSG:4326"){
  print("warping raster...")
  srcfile <- gsub("[\\/]$", "", srcfile)

  if (missing(dstfile)) dstfile <- gsub("\\.(.{2,5})$", "_warped\\.sdat", srcfile)

  gdalUtils::gdalwarp(srcfile = srcfile,
           dstfile = dstfile,
           srcnodata = srcnodata, dstnodata = dstnodata,
           s_srs = s_srs, t_srs = output_crs, of = "SAGA", r = "near",
           overwrite = T,  stderr = T)

  return(dstfile)
}

#===============================================================================
#' @title Mosaic rasters with GDAL
#===============================================================================
gdal_mosaic <- function(srcfiles, dstfile, srcnodata = -32768, dstnodata, of = "SAGA"){
  print("building mosaic...")
  gdalUtils::mosaic_rasters(gdalfile = srcfiles,
                      dst_dataset = dstfile,
                      srcnodata = srcnodata, dstnodata = dstnodata,
                      of = "SAGA", r = "near", overwrite = T, stderr = T)

  return(dstfile)
}

#===============================================================================
#' @title Mosaic and warp raster set
#'
#' @description combines some grids and changes coordinate system
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @param output_crs character vector specifying EPSG for coordinate system to
#' output
#'
#' @keywords internal
#'
#' @export
#===============================================================================
MosaicAndWarp <- function(gridnames, DEM_path, saga_env,
                  input_crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                  output_crs="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96
                  +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"){
  #TODO: this could be made faster I expect (maybe use the point as mosaic bounds?)
  # First, convert tiles to sgrd format so that saga can use them
  # (keep original CRS),  then mosaic the tiles using SAGA,
  # and finally, warp the mosaicked grid using gdal into the desired CRS

  temp.files <- character()

  for (i in seq(1, length(gridnames))){
    x <- gridnames[i]
    outputfile <- file.path(saga_env$workspace,
                            paste("temp", i, ".sdat", sep = "")) # make dummy temp file

    input.DEM <- GetTilePathsHS(x, DEM_path)
    temp.files <- c(temp.files, outputfile)
    gdal_warp2SAGA(input.DEM, outputfile, output_crs = input_crs)
  }

  temp_hdrs_path <- gsub("sdat", "sgrd", temp.files)
  temp_hdrs <- lapply(temp_hdrs_path, read.sgrd.header)
  temp_hdrs <- do.call(rbind, temp_hdrs)

  minx <- as.character(floor(min(as.numeric(temp_hdrs$POSITION_XMIN))))
  miny <- as.character(floor(min(as.numeric(temp_hdrs$POSITION_YMIN))))
  maxx <- as.character(max(as.numeric(temp_hdrs$POSITION_XMIN) +
             as.numeric(temp_hdrs$CELLCOUNT_X) * as.numeric(temp_hdrs$CELLSIZE)))

  maxy <- as.character(max(as.numeric(temp_hdrs$POSITION_YMIN) +
             as.numeric(temp_hdrs$CELLCOUNT_Y) * as.numeric(temp_hdrs$CELLSIZE)))

  cellsize <- temp_hdrs$CELLSIZE[1]


  #mosaic
  inputs <- paste(temp_hdrs_path, collapse = ";")
  mosaicked <- file.path(saga_env$workspace, "temp_mosaicked.sdat")
  MoisaicRS(inputs, mosaicked, saga_env = saga_env, xmin = minx, xmax = maxx,
            ymin = miny, ymax = maxy, cellsize = cellsize)

  warped.output <-  file.path(saga_env$workspace, "temp_mosaicked_warped.sdat")
  gdal_warp2SAGA(mosaicked, warped.output, output_crs = output_crs)

  return(warped.output)
}

#===============================================================================
#' @title Build channel network using RSAGA
#'
#' @description Builds a channel network from a DEM.
#'
#' @param in_DEM character string pointing to the *.sgrd path of a
#' (hydrologically appropriate) DEM
#'
#' @param init.grid (optional) character string pointing to a flow accumulation
#' *.sgrd file.  If not supplied, it will be created.
#'
#' @param out.shape character string naming a shapefile open for writing which
#' will be the channel network
#'
#' @param upstr_threshold integer, number of upstream cells required to
#' initiate channel growth. A smaller number will yield more channels and will
#' snap to smaller streams. Larger numbers will produce only major channel
#'  branches. For a 50k DEM, numbers between 20,000 and 50,000 are appropriate
#'  to reproduce the main channel network without many small branches.
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' character string  Mosaic 'grid_tools' module 3
#'
#' @return character string path to channel network shapefile
#'
#' @export
#===============================================================================
ChannelNetworkRS <- function(in_DEM, init.grid, out.shape,
                             upstr_threshold = 35000, saga_env, verbose = F){


  if (missing(init.grid)){
    init.grid <- FlowAccumulationRS(in_DEM, saga_env = saga_env)
  }

  if (missing(out.shape)){
    out.shape <- gsub("\\.s[gd][ra][dt]", "_channels\\.shp", in_DEM)
  }

  print("Building Channel Network...")

  RSAGA::rsaga.geoprocessor(lib = "ta_channels", module = 0, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
                              ELEVATION   = in_DEM,
                              INIT_GRID   = init.grid,
                              SHAPES      = out.shape,
                              INIT_METHOD = 2,
                              INIT_VALUE  = upstr_threshold
                            ))
  return(out.shape)
}

#===============================================================================
#' @title Generate Flow Accumulation using RSAGA
#'
#' @description Calculates flow accumulation
#' @param in_DEM character string pointing to the *.sgrd path of a
#' (hydrologically appropriate) DEM
#' @param output_grid character string pointing to the *.sgrd path of the
#' output file
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @return path to output flow accumulation grid
#' @export
#===============================================================================
FlowAccumulationRS <- function(in_DEM, out_grid, saga_env, verbose = F){
  if (missing(out_grid)) out_grid <- gsub("\\.sgrd", "_flowacc\\.sgrd", in_DEM)
  print("Calculating Flow Accumulation...")

  RSAGA::rsaga.geoprocessor(lib = "ta_hydrology", module = 0, env = saga_env,
                            show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
                      ELEVATION  = in_DEM,
                      CAREA      = out_grid,  # listed as FLOW and FLOW_UNIT in the GUI
                      CAREA_UNIT = 0, # number of cells (less accurate but faster)
                      METHOD     = 0     # deterministic 8
                            ))

  return(out_grid)
}

#===============================================================================
#' @title Measure raster at point using RSAGA
#'
#' @description Spatial Join: Retrieves information from the selected grid at
#' the positions of the points and returns the grid values.
#'
#' @param point either SpatialPoints* or character path to points used to
#' intersect grid. Must have same CRS as the grid
#'
#' @param grid grid whose values will be sampled
#'
#' @param saga_env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#'
#' @return vector of grid values the same length as the number of points
#'
#' @export
#===============================================================================
SampleRasterRS <- function(point, grid, saga_env, verbose=F){
  print("Sampling raster...")

  if (grepl("Spatial", class(point))){
    rgdal::writeOGR(point, dsn = file.path(saga_env$workspace, "station.shp"),
                    layer = "station.shp", driver = "ESRI Shapefile")
    point <- file.path(saga_env$workspace, "station.shp")
  }

  sampled <- gsub("\\.shp", "_gridsample\\.shp", point)

  RSAGA::rsaga.geoprocessor(lib = "shapes_grid", module = 0, env = saga_env,
                            show.output.on.console = verbose,
                     check.module.exists = FALSE,
                     param = list(
                       SHAPES   = point,
                       GRIDS    = grid,
                       RESULT   = sampled,
                       INTERPOL = 0
                     ))

  x <- invisible(rgdal::readOGR(sampled))
  x <- x@data[, ncol(x@data)]

  return(x)
}

#===============================================================================
#' @title Calculate upslope polygon from DEM
#'
#' @description Create a drainage basin shapefile using a DEM and a 'pour point'.
#'  The polygon is
#'
#' @param point a SpatialPointsDataFrame corresponding to a hydrometric station.
#' Must be in same coordinates system as DEM and coordinate system must be
#' projected (not lat/long).  Must have longitude and latitude attributes.
#'
#' @param saga_env an rsaga environment object
#'
#' @param outdir Directory to output final upslope shapefile
#'
#' @param pointbuffer numeric, how much should the point be buffered (calculates
#' upslope basin of buffered area)
#'
#' @param iterate_to numeric, if greater than pointbuffer, checks that the
#' calculated upslope area is at least 0.5 sq. km larger than the
#'
#' @param iterate_incr numeric, by how much should the buffer be grown each
#' iteration (ignored if iterate_to is equal to pointbuffer)
#'
#'  @param iterate_thres numeric, how many square kilometers difference between
#'  the calculated area and initial point buffer should be required to trigger
#'  another iteration (ignored if iterate_to is equal to pointbuffer)
#'
#' @param DEM_source character, one of: c('CDED', 'NED', 'CDEM, CDSM', 'SHEDS').
#' Ignored if a DEM file is supplied to DEM_path
#'
#' @param DEM_path One of:
#'  (1) a file path to a directory containing DEM files either in the the format
#'  n%%w0%%_con_grid.sgrd (if DEM_source is 'SHEDS') or a directory to which DEM
#'  tiles will be downloaded (if DEM_source='CDEM', 'CDED' etc.),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a
#'  projected coordinate system and the coordinatesystem should match that of
#'  the point (e.g. Canada Albers Conformal Conic)
#'
#'  @param clip_size How big (in m) a clip should be generated from the
#'  original DEM (too big doesn't work and is slower)
#'
#'  @param projected.CRS character vector (proj4) specifying "working" CRS in
#'  which area and distance calculations are to be done
#'
#' @param ... optional arguments to \link[RSAGA]{rsaga.geoprocessor} including
#' 'verbose' (toggles SAGA console output)
#'
#' @param method character, one of ('D8', 'DINF'), specifying the hydrological
#' pathing model
#'
#' @param upstr_threshold integer, number of upstream cells required to
#' initiate channel growth. A smaller number will yield more channels and will
#'  snap to smaller streams. Larger numbers will produce only major channel
#'   branches. For a 50k DEM, numbers between 20,000 and 50,000 are appropriate
#'    to reproduce the main channel network without many small branches.
#'
#'  @return A list containing diagnostic information and and a character string
#'   naming the recently created upstream area polygon
#'
#' @export
#===============================================================================
UpslopeDEM <- function(point, DEM_path, DEM_source = "NTS", saga_env, outdir,
                       pointbuffer = 50, clip_size = 50000,
                       projected.CRS, iterate_to = 150, iterate_incr = 50,
                       iterate_thres = 0.4, method = "D8",
                       upstr_threshold = 35000, ...){

  # set some parameters
  method  <- switch(toupper(method), "D8" = 0, "DINF" = 1)
  snapped <- F
  pb_init <- pointbuffer

  #sanity checks
  if (missing(projected.CRS)){
    projected.CRS <- GetProj4("AlbersEqualAreaConic")
  }
  iterate_to <- max(c(iterate_to, pointbuffer))
  in_DEM <- gsub("\\.sdat$", "\\.sgrd", DEM_path)

  # prepare workspace
  oldwd <- getwd()
  setwd(saga_env$workspace)
  outdir <- gsub("[\\/]$", "", outdir)

  # get name of station
  station_name <- point@data[1, grepl("^\\d{2}[[:alpha:]]{2}\\d{3}$", point@data)]
  print(station_name)
  final_name <- file.path(outdir, paste(station_name, "_upslope.shp", sep = ""))

  #if point not in proper crs, convert it
  if (point@proj4string@projargs != projected.CRS){
    point <- sp::spTransform(point, CRSobj = sp::CRS(projected.CRS))
  }

  # Find missing DEM if necessary / Convert to sgrd if necessary.
  if (!grepl("\\.sgrd$", in_DEM)){
    # if not an SGRD file (then we expect a directory or a list of names)

    if (tolower(DEM_source) == "sheds"){
      name <- HydroMosaic(point@data$longitude, point@data$latitude,
                          tol = clip_size)

      if (length(name) == 1){
        matching_files <- list.files(in_DEM,
                                      pattern = paste(name, ".sgrd", sep = ""),
                                      full.names = T, recursive = T)

        if (length(matching_files) == 0){
          #  if sgrid doesn't exist, create it
          original.DEM <- GetTilePathsHS(name, in_DEM)
          print(sprintf("Converting %s DEM to sgrd...", original.DEM))
          gdal_warp2SAGA(original.DEM, output_crs = projected.CRS)
        }

        in_DEM <- list.files(in_DEM, pattern = paste(name, ".sgrd$", sep = ""),
                             full.names = T, recursive = T)

        print(sprintf("Using %s as input DEM", in_DEM))

      }else if (length(name) > 1){
        # need to mosaic and transform
        print("Multiple grids Necessary...")
        original.DEM <- GetTilePathsHS(name, in_DEM)
        in_DEM <- MosaicAndWarp(gridnames = name, DEM_path = in_DEM,
                                saga_env = saga_env, output_crs = projected.CRS)
      }

    }else if (toupper(DEM_source) %in% c("CDED", "NED", "CDEM", "CDSM")){
      print(DEM_source) #debug only
      in_DEM <- OverlayDEM(point, DEM_dir = DEM_path,
                           output.dir = saga_env$workspace,
                           product = DEM_source,
                           tol = clip_size) # clip from DEM

      clip_size <- 0
    }else{
      print("Could not find suitable DEM in folder")

      return()
    }
  }

  # Test if point is outside of raster limits (e.g. in USA)
  value <- SampleRasterRS(point = point, grid = in_DEM, saga_env = saga_env )

  nodata <- F
  if (is.na(value)){
    print("station outside of CDED coverage, switching to NED data")
    in_DEM <- OverlayDEM(point, DEM_dir = DEM_path,
                         output.dir = saga_env$workspace, product = "NED",
                         tol = clip_size)
    # nodata <- T
    # setwd(oldwd)
    # return(list(final_name=NA, pointbuffer=NA, iterate_thres=iterate_thres,
    #             iterate_to=iterate_to,
    #             nodata=nodata, snapped=F, snapped=snapped, snap_dist=NA))
  }

  # Clip grid to point
  if (clip_size > 0){
    ClipGridRS(point, in_DEM, "clipped.sgrd", saga_env,
               tol = clip_size, ...)
  }

  # Fill Sinks
  FillSinksRS(ifelse( (clip_size > 0), "clipped.sgrd", in_DEM),
              "filled.sgrd", saga_env, MINSLOPE = 0.01, ...)

  snapped <- FALSE
  snap_dist <- NA
  run <- TRUE

  while (run){
    # Point to Grid
    Point2GridRS(point, "filled.sgrd", "target.sgrd",  saga_env,
                 buffer_dist = pointbuffer, ...)

    # Upslope
    upslope <- UpslopeAreaRS("filled.sgrd", "target.sgrd", "upslope.sgrd",
                             saga_env, method = method, ...)
    if (method == 1){
      upslope <- GridThresholdRS(upslope,
                                 threshold_value = 0,
                                 saga_env = saga_env)
    }

    # Find Area
    area <- GridVolumeRS(upslope, level = 99, saga_env = saga_env) #  (in m2)

    if (abs(area - pi * pointbuffer ** 2) * 1e-6 < iterate_thres |
        pointbuffer > iterate_to ){

      # if upslope area didn't work, try widening search radius
      pointbuffer <- pointbuffer + iterate_incr

      if (pointbuffer > iterate_to | snapped == TRUE){
        if (snapped == TRUE){
          run <- FALSE
        }else{
          snapped <- TRUE
          snapped_pt <- SnapToPourPoint(point, "filled.sgrd",
                                        saga_env = saga_env,
                                        upstr_threshold = upstr_threshold)
          point <- snapped_pt$station
          snap_dist <- snapped_pt$distance
          pointbuffer <- pb_init
        }
      }
    }else{
      run <- FALSE
    }
  }

  # Convert to polygon
  Grid2PolyRS(upslope, final_name, saga_env, ...)

  # put workspace back the way it was
  setwd(oldwd)
  return(list(final_name = final_name, pointbuffer = pointbuffer,
              iterate_thres = iterate_thres, iterate_to = iterate_to,
              nodata = nodata, snapped = snapped, snap_dist = snap_dist))
}
