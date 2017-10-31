#' troubleshooting RSAGA:
#' - run command in commandline
#' - sometimes the variables names are different between versions of SAGA, (e.g. INTERPOL instead of RESAMPLING)




#' Upslope polygon
#'
#' @description Create a drainage basin shapefile using a DEM and a 'pour point'.  The polygon is
#' @param point a SpatialPointsDataFrame corresponding to a hydrometric station.  Must be in same coordinates system as DEM and coordinate system
#' must be projected (not lat/long).  Must have longitude and latitude attributes.
#' @param saga.env an rsaga environment object
#' @param outdir Directory to output final upslope shapefile
#' @param pointbuffer numeric, how much should the point be buffered (calculates upslope basin of buffered area)
#' @param iterate.to numeric, if greater than pointbuffer, checks that the calculated upslope area is at least 0.5 sq. km larger
#' than the
#' @param iterate.incr numeric, by how much should the buffer be grown each iteration (ignored if iterate.to is equal to pointbuffer)
#' @param iterate.thres numeric, how many square kilometers difference between the calculated area and
#' initial point buffer should be required to trigger another iteration (ignored if iterate.to is equal to pointbuffer)
#' @param DEM.path One of:
#'  (1) a file path to a directory containing DEM files either in the the format n%%w0%%_con_grid.sgrd (if DEM.source is 'SHEDS')
#'  or a directory to which NTS files will be downloaded (if DEM.source='NTS'),
#'  (2) a file path to a dem file in SAGA format.  The DEM should be in a projected coordinate system and the coordinate
#'  system should match that of the point (e.g. Canada Albers Conformal Conic)
#' @param dem.clip.square How big (in m) a clip should be generated from the original DEM (too big doesn't work and is slower)
#' @param projected.CRS character vector (proj4) specifying "working" CRS in which area and distance calculations are to be done
#' @param ... optional arguments including 'verbose' (toggles SAGA console output)
#' @param method character, one of ('D8', 'DINF'), specifying the hydrological pathing model
#' @param upstream.threshold integer, number of upstream cells required to initiate channel growth.
#' A smaller number will yield more channels and will snap to smaller streams. Larger numbers will produce
#' only major channel branches. For a 50k DEM, numbers between 20,000 and 50,000 are appropriate to reproduce the
#' main channel network without many small branches.
#' @return Character string naming the recently created polygon
#' @export
UpslopeDEM <- function(point, DEM.path, DEM.source='NTS', saga.env, outdir, pointbuffer=50, dem.clip.square=50000,
                       projected.CRS, iterate.to=150, iterate.incr=50, iterate.thres=0.4, method="D8", upstream.threshold=35000, ...){

  # set some parameters
  method <- switch(toupper(method), "D8"=0, "DINF"=1)
  snapped<-F
  pb.init <- pointbuffer

  #sanity checks
  if (missing(projected.CRS)){
    projected.CRS <- GetProj4("AlbersEqualAreaConic")
  }
  iterate.to <- max(c(iterate.to, pointbuffer))
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

    if (grepl('sheds',tolower(DEM.source))){
      name <- HydroMosaic(point@data$longitude, point@data$latitude, tol=dem.clip.square)
      if (length(name)==1){
        check.if.exists <- list.files(in.DEM, pattern=paste(name, ".sgrd", sep=''), full.names = T, recursive = T)
        if (length(check.if.exists) ==0){  #  if sgrid doesn't exist, create it
          original.DEM <- GetTilePathsHS(name, in.DEM)
          print(sprintf("Converting %s DEM to sgrd...", original.DEM))
          gdal_warp2SAGA(original.DEM, outputCRS = projected.CRS)
        }
        in.DEM <- list.files(in.DEM, pattern=paste(name,".sgrd$", sep=''), full.names = T, recursive = T)
        print(sprintf("Using %s as input DEM", in.DEM))

      }else if (length(name) > 1){ # need to mosaic and transform
        print('Multiple grids Necessary...')
        original.DEM <- GetTilePathsHS(name, in.DEM)
        in.DEM <- MosaicAndWarp(gridnames = name, DEM.path = in.DEM, saga.env = saga.env, outputCRS = projected.CRS)
      }

      }else if (grepl('nts',tolower(DEM.source))){
      in.DEM <- OverlayNTS(point, NTS.dir=DEM.path, output.dir=saga.env$workspace, tol=dem.clip.square) # clip from DEM
      dem.clip.square <- 0
    }else{
      print("Could not find suitable DEM in folder")
      return()
    }
  }

  # If point is outside of raster limits (e.g. in USA)
  value <- SampleRasterRS(point = point, grid = in.DEM, saga.env=saga.env )
  nodata <- F
  if (is.na(value)){
   nodata <- T
   setwd(oldwd)
   return(list(final.name=NA, pointbuffer=NA, iterate.thres=iterate.thres, iterate.to=iterate.to,
               nodata=nodata, snapped=F, snapped=snapped, snap.dist=NA))
  }

  # Clip grid to point
  if (dem.clip.square > 0){
  ClipGridRS(point, in.DEM, 'clipped.sgrd', saga.env, tol=dem.clip.square, ...)
  }

  # Fill Sinks
  FillSinksRS(ifelse((dem.clip.square > 0),'clipped.sgrd',in.DEM),
              'filled.sgrd', saga.env, MINSLOPE=0.01, ...)

snapped <- FALSE
snap.dist <- NA
run <- TRUE
  while (run){
    # Point to Grid
    Point2GridRS(point, 'filled.sgrd', 'target.sgrd',  saga.env, buffer.dist=pointbuffer, ...)

    # Upslope
    upslope <- UpslopeAreaRS('filled.sgrd', 'target.sgrd', 'upslope.sgrd', saga.env, method=method, ...)
    if (method==1){
      upslope <- GridThresholdRS(upslope, threshold.value = 0, saga.env=envi)
    }

    # Find Area
    area <- GridVolumeRS(upslope, level = 99, saga.env=envi) # returns in m2

    if (abs(area - pi*pointbuffer**2)*1e-6 < iterate.thres){  # if upslope area didn't work, try widening search radius
      pointbuffer <- pointbuffer + iterate.incr

    if (pointbuffer > iterate.to | snapped==TRUE){
      if (snapped==TRUE){
        run <- FALSE
      }else{
        snapped <- TRUE
        snapped.pt <- SnapToPourPoint(point, 'filled.sgrd', saga.env=saga.env, upstream.threshold=upstream.threshold)
        point <- snapped.pt$station
        snap.dist <- snapped.pt$distance
        pointbuffer <- pb.init
      }
    }
    }else{
      run <- FALSE
    }
  }

  # Convert to polygon
  Grid2PolyRS(upslope, final.name, saga.env, ...)

  # put workspace back the way it was
  setwd(oldwd)
  return(list(final.name=final.name, pointbuffer=pointbuffer, iterate.thres=iterate.thres, iterate.to=iterate.to,
              nodata=nodata, snapped=snapped, snap.dist=snap.dist))
}


#' @description
#' @param station spatialpointsdataframe representing a hydrological station
#' @param in.DEM character path pointing to hydrologically conditioned DEM. Must be in projected
#' coordinate system and, at a minimum, have filled sinks.
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param upstream.threshold integer, number of upstream cells required to initiate channel growth.
#' A smaller number will yield more channels and will snap to smaller streams. Larger numbers will produce
#' only major channel branches. For a 50k DEM, numbers between 20,000 and 50,000 are appropriate to reproduce the
#' main channel network without many small branches.
#' @return a SpatialPointsDataFrame with the same data as the input station, but located on a DEM stream branch.
#' @export
SnapToPourPoint <- function(station, in.DEM, saga.env, upstream.threshold=35000){
  chann <- ChannelNetworkRS(in.DEM=in.DEM, upstream.threshold=upstream.threshold, saga.env=envi, verbose=F)
  chann <- rgdal::readOGR(chann)
  newpoint <- SnapToNearest(station, chann)
  dist <- rgeos::gDistance(station, chann)
  return(list(station=newpoint, distance=dist))
}


#' Shapes to Grid
#'
#' @description Turns a vector object into a grid ('grid_gridding' module 0)
#' @param point A SpatialPoint or SpatialPointsDataFrame corresponding to a station
#' @param grid.sys A SAGA *.sgrd file with the same grid system as the desired output
#' @param out.name Character string of output file name
#' @param env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param buffer.dist If equal to zero does not buffer.
#' @keywords internal
#' @export
Point2GridRS <- function(point, grid.sys, out.name,  saga.env, buffer.dist=150, verbose=F){
  print("Converting point to grid ...")
  # Write point as shapefile
  rgdal::writeOGR(point, dsn=file.path(saga.env$workspace, "target.shp"), layer="target.shp", driver="ESRI Shapefile")
  # Convert shapefile to grid
  RSAGA::rsaga.geoprocessor(lib = 'grid_gridding', module=0, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    INPUT=file.path(saga.env$workspace, "target.shp"),
    OUTPUT=0,
    GRID_TYPE=0,
    TARGET_DEFINITION=1,  # Better to define by the other grid (less prone to error? but maybe slower??)
    TARGET_TEMPLATE=grid.sys,
    GRID=out.name
  ))

  if (buffer.dist > 0){
    # Grid buffer module 8
    RSAGA::rsaga.geoprocessor(lib = 'grid_tools', module=8, env=saga.env, show.output.on.console = verbose,
                              check.module.exists = FALSE, param = list(
      FEATURES=out.name,
      BUFFER=out.name,
      DIST=buffer.dist,
      BUFFERTYPE=0
    ))
  }
}

#' Clip grid to point plus error  # Clip Grids = 31
#'
#' @param point SpatialPointsDataFrame
#' @param in.DEM character path to tiled grid to clip
#' @param out.DEM (optional) path to save output DEM, otherwise appends '_clipped' to filename of in.DEM
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param tol size of clipping box (metres) around point
#' @keywords internal
#' @export
ClipGridRS <- function(point, in.DEM, out.DEM, saga.env, tol, verbose=F){
  print("Clipping grid based on point location...")
  if (missing(out.DEM)) out.DEM <- gsub("\\.sgrd", "_clipped\\.sgrd", in.DEM)
  x = point@coords[,1]
  y = point@coords[,2]

  RSAGA::rsaga.geoprocessor(lib = 'grid_tools', module=31, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    GRIDS=in.DEM,
    CLIPPED=out.DEM,
    XMIN=x - tol,
    XMAX=x + tol,
    YMIN=y - tol,
    YMAX=y + tol
  ))
}

#' Clip Polygon   # Clip Grids = 11
#'
#' @param input SpatialPointsDataFrame
#' @param clipping.layer character path to tiled grid to clip
#' @param output (optional) path to save output DEM, otherwise appends '_clipped' to filename of in.DEM
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @keywords internal
#' @export
ClipPolygonRS <- function(input, clipping.layer, output, saga.env,verbose=F){
  print("Clipping Polygons...")
  if (missing(output)) output <- gsub("\\.sgrd", "_clipped\\.sgrd", input)

  RSAGA::rsaga.geoprocessor(lib = 'shapes_polygons', module=11, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
                              CLIP=clipping.layer,
                              S_INPUT=input,
                              S_OUTPUT=output,
                              DISSOLVE=1
                              ))
}

#' Fill Sinks
#'
#' @description Fill sinks according to Wang and Liu ('ta_preprocessor' module 4 or 5 for XXL)
#' @param in.DEM character path to grid to fill
#' @param out.DEM (optional) path to save output DEM, otherwise appends '_filled' to filename of in.DEM
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param MINSLOPE minimum slope to maintain between cells
#' @param XXL logical, whether or not to use the 'Fill sinks XXL' tool, suitable for large grids
#' @keywords internal
#' @export
FillSinksRS <- function(in.DEM, out.DEM, saga.env, MINSLOPE=0.01, verbose=F, XXL=F){
  print("Filling sinks...")
  if (missing(out.DEM)) out.DEM <- gsub("\\.s[gd][ra][dt]", "_filled\\.sgrd", in.DEM)

  # Fill sinks (wang & Liu) = 4
  module <- ifelse(XXL, 5, 4)
  RSAGA::rsaga.geoprocessor(lib = 'ta_preprocessor', module=module, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param=list(
                              ELEV=in.DEM,
                              FILLED=out.DEM,
                              MINSLOPE=MINSLOPE  # not sure how to choose optimal value
                            ))
  return(out.DEM)
}

#' Sink Removal
#'
#' @description Remove sinks by filling sinks or deepening drainage routes using the method of Conrad (2001)
#' @param in.DEM character path to grid to fill
#' @param out.DEM (optional) path to save output DEM, otherwise appends '_filled' to filename of in.DEM
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param METHOD Available Choices: [0] Deepen Drainage Routes [1] Fill Sinks.
#' @param THRESHOLD logical, whether or not to exclude deeper sinks from filling
#' @param THRESHOLD.HEIGHT numeric, The parameter describes the maximum depth of a sink to be considered for removal [map units].
#'  This makes it possible to exclude deeper sinks from filling.
#' @keywords internal
#' @export
BreachDepressionsRS <- function(in.DEM, out.DEM, saga.env, METHOD=0, THRESHOLD=0, THRESHOLD.HEIGHT=100, verbose=F){
  print(paste(ifelse(METHOD, "Filling", "Breaching"), "sinks..."))

  if (missing(out.DEM)) out.DEM <- gsub("\\.sgrd", "_sinksremoved\\.sgrd", in.DEM)


  RSAGA::rsaga.geoprocessor(lib = 'ta_preprocessor', module=2, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param=list(
                              DEM=in.DEM,
                              DEM_PREPROC=out.DEM,
                              METHOD=METHOD,
                              THRESHOLD=THRESHOLD,
                              THRSHEIGHT=THRESHOLD.HEIGHT
                            ))
  return(out.DEM)
}

#' Upslope area
#'
#' @description Calculates upslope area and produces a grid ('ta_hydrology' module 4)
#' @param in.DEM Hydrologically appropriate DEM file (saga grid)
#' @param target.grid saga grid
#' @param out.GRD character string path to output grid
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @keywords internal
#' @export
UpslopeAreaRS <- function(in.DEM, target.grid, out.GRD, saga.env, verbose=F, method=1){
  print("Calculating upstream area...")
  RSAGA::rsaga.geoprocessor(lib='ta_hydrology', module=4, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param=list(
    TARGET=target.grid,
    ELEVATION=in.DEM,
    AREA=out.GRD,
    METHOD=method
  ))
  return(out.GRD)
}

#' Grid to polygon
#'
#' @description Creates a polygon from a grid ('shapes_grid' module 6: Vectorising Grid Classes)
#' @param in.GRD character string path to SAGA grid file.  The grid should
#' @param out.poly character string path to output polygon shapefile
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param CLASS_ID numeric raster value that designates which grid cells to convert to polygon (cells with other values will be ignored)
#' @keywords internal
#' @export
Grid2PolyRS <- function(in.GRD, out.poly, saga.env, CLASS_ID=100, verbose=F){
  print("Converting grid to polygons...")
  RSAGA::rsaga.geoprocessor(lib = 'shapes_grid', module=6, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    GRID=in.GRD,
    POLYGONS=out.poly,
    CLASS_ALL=0,
    CLASS_ID=CLASS_ID,
    SPLIT=1,
    ALLVERTICES=0
  ))
}

# 'shapes_polygons' module 15
DifferenceRS <- function(shape, eraser, output, saga.env, verbose=F){
  print("Differencing shapes...")
  RSAGA::rsaga.geoprocessor(lib = 'shapes_polygons', module=15, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    A=shape,
    B=eraser,
    RESULT=output,
    SPLIT=0
    ))
}

# 'grid_calculus' module 2
GridVolumeRS <- function(grid, level, method=0, saga.env){
  print("Calculating Grid Volume...")
  volume <- invisible(capture.output(RSAGA::rsaga.geoprocessor(lib = 'grid_calculus', module=2, env=saga.env,
                                                               check.module.exists = FALSE, param = list(
                              GRID=grid,
                              METHOD=method,
                              LEVEL=level
                            ))))
  volume <- as.numeric(gsub("Volume: ", "", regmatches(volume,regexpr("^Volume: (\\d*\\.\\d*)", volume))))
  return(volume)
}

GridThresholdRS <- function(grid, threshold.value, method=0, saga.env, dstfile){
  print("Thresholding Grid...")
  if (missing(dstfile)) dstfile <- gsub("\\.(.{2,5})$","_thres\\.sdat", grid)
  volume <- capture.output(RSAGA::rsaga.geoprocessor(lib = 'grid_calculus', module=1, env=saga.env,
                                                     check.module.exists = FALSE, param = list(
                                                                 GRIDS=grid,
                                                                 FORMULA=sprintf("(g1>%d)*100",threshold.value),
                                                                 TYPE=4,
                                                                 RESULT=dstfile
                                                                  )))

  return(dstfile)
}

# Mosaic 'grid_tools' module 3
MoisaicRS <- function(grids, out_grid, saga.env, xmin, xmax, ymin, ymax, cellsize, verbose=F){
  print("Building Mosaic...")
  RSAGA::rsaga.geoprocessor(lib = 'grid_tools', module=3, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
    GRIDS=grids,
    TYPE=4,  # 2-byte signed integer
    INTERPOL=4,
    OVERLAP=0,
    MATCH=0,
    TARGET_DEFINITION=0,
    TARGET_USER_XMIN=xmin,
    TARGET_USER_XMAX=xmax,
    TARGET_USER_YMIN=ymin,
    TARGET_USER_YMAX=ymax,
    TARGET_USER_SIZE=cellsize,
    TARGET_OUT_GRID=out_grid
    ))
}



#' Convert to SAGA grid
#'
#' @description Convertes a GRID DEM into SAGA .sgrd/.sdat format.  Also converts the CRS into a projected
#' system.
#' @param srcfile path to GRID folder (no trailing path separators)
#' @param dstfile (optional) path to output SAGA grid file.  If not specified, saves to same directory as input.  Note
#' that the file name must have the "sdat" file extension (not sgrid)
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param outputCRS character string specifying the coordinate system of the output grid.  Examples include:
#' "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
#'
#' @keywords internal
#' @export
gdal_warp2SAGA <- function(srcfile, dstfile, outputCRS, srcnodata=-32768, dstnodata=-32767, s_srs="EPSG:4326"){
  print("warping raster...")
  srcfile <- gsub('[\\/]$', '', srcfile)
  if (missing(dstfile)) dstfile <- gsub("\\.(.{2,5})$","_warped\\.sdat",srcfile)
  gdalUtils::gdalwarp(srcfile=srcfile,
           dstfile=dstfile,
           srcnodata=srcnodata, dstnodata=dstnodata,
           s_srs = s_srs, t_srs=outputCRS, of = "SAGA", r='near', overwrite=T, stderr=T)
  return(dstfile)
}

gdal_mosaic <- function(srcfiles, dstfile, srcnodata=-32768, dstnodata, of="SAGA"){
  print("building mosaic...")
  gdalUtils::mosaic_rasters(gdalfile=srcfiles,
                      dst_dataset=dstfile,
                      srcnodata=srcnodata, dstnodata=dstnodata,
                      of = "SAGA", r='near', overwrite=T, stderr=T)
  return(dstfile)
}

#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' @param outputCRS character vector specifying EPSG for coordinate system to output
#' @keywords internal
#' @export
MosaicAndWarp <- function(gridnames, DEM.path, saga.env, inputCRS="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs",
                          outputCRS="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"){  #TODO: this could be made faster I expect (maybe use the point as mosaic bounds?)
  # First, convert tiles to sgrd format so that saga can use them (keep original CRS), then mosaic the tiles using SAGA,
  # and finally, warp the mosaicked grid using gdal into the desired CRS

  temp.files <- character()
    for (i in seq(1,length(gridnames))){
    x <- gridnames[i]
    outputfile <- file.path(saga.env$workspace, paste("temp",i,'.sdat',sep='')) # make dummy temp file
    input.DEM <- GetTilePathsHS(x, DEM.path)
    temp.files <- c(temp.files, outputfile)
    gdal_warp2SAGA(input.DEM, outputfile, outputCRS = inputCRS)
  }
  temp.hdrs.path <- gsub("sdat", "sgrd", temp.files)
  temp.hdrs <- lapply(temp.hdrs.path, read.sgrd.header)
  temp.hdrs <- do.call(rbind, temp.hdrs)

  minx <- as.character(floor(min(as.numeric(temp.hdrs$POSITION_XMIN))))
  miny <- as.character(floor(min(as.numeric(temp.hdrs$POSITION_YMIN))))
  maxx <- as.character(max(as.numeric(temp.hdrs$POSITION_XMIN) + as.numeric(temp.hdrs$CELLCOUNT_X) * as.numeric(temp.hdrs$CELLSIZE)))
  maxy <- as.character(max(as.numeric(temp.hdrs$POSITION_YMIN) + as.numeric(temp.hdrs$CELLCOUNT_Y) * as.numeric(temp.hdrs$CELLSIZE)))
  cellsize <- temp.hdrs$CELLSIZE[1]


  #mosaic
  inputs <- paste(temp.hdrs.path, collapse=';')
  mosaicked <- file.path(saga.env$workspace, "temp_mosaicked.sdat")
  MoisaicRS(inputs, mosaicked, saga.env=saga.env, xmin=minx, xmax=maxx, ymin=miny, ymax=maxy, cellsize=cellsize)
warped.output <-  file.path(saga.env$workspace, "temp_mosaicked_warped.sdat")
  gdal_warp2SAGA(mosaicked, warped.output, outputCRS = outputCRS)
  return(warped.output)
}

#' @description Builds a channel network.
#' @param in.DEM character string pointing to the *.sgrd path of a (hydrologically appropriate) DEM
#' @param init.grid (optional) character string pointing to a flow accumulation *.sgrd file.  If not
#' supplied, it will be created.
#' @param out.shape character string naming a shapefile open for writing which will be the channel network
#' @param upstream.threshold integer, number of upstream cells required to initiate channel growth.
#' A smaller number will yield more channels and will snap to smaller streams. Larger numbers will produce
#' only major channel branches. For a 50k DEM, numbers between 20,000 and 50,000 are appropriate to reproduce the
#' main channel network without many small branches.
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
#' character string  Mosaic 'grid_tools' module 3
ChannelNetworkRS <- function(in.DEM, init.grid, out.shape, upstream.threshold=35000, saga.env, verbose=F){
  print("Building Channel Network...")
  if (missing(init.grid)){
    init.grid <- FlowAccumulationRS(in.DEM, saga.env=saga.env)
  }
  if (missing(out.shape)){
    out.shape <- gsub("\\.s[gd][ra][dt]", "_channels\\.shp", in.DEM)
  }
  RSAGA::rsaga.geoprocessor(lib = 'ta_channels', module=0, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
                              ELEVATION=in.DEM,
                              INIT_GRID=init.grid,
                              SHAPES=out.shape,
                              INIT_METHOD=2,
                              INIT_VALUE=upstream.threshold
                            ))
  return(out.shape)
}

#' @description Calculates flow accumulation
#' @param in.DEM character string pointing to the *.sgrd path of a (hydrologically appropriate) DEM
#' @param output_grid character string pointing to the *.sgrd path of the output file
#' @param saga.env A SAGA geoprocessing object.  Suggested version is 2.2.2.
FlowAccumulationRS <- function(in.DEM, out_grid, saga.env, verbose=F){
  if (missing(out_grid)) out_grid <- gsub("\\.sgrd", "_flowacc\\.sgrd", in.DEM)
  print("Calculating Flow Accumulation...")
  RSAGA::rsaga.geoprocessor(lib = 'ta_hydrology', module=0, env=saga.env, show.output.on.console = verbose,
                            check.module.exists = FALSE, param = list(
                              ELEVATION=in.DEM,
                              CAREA=out_grid,  # listed as FLOW and FLOW_UNIT in the GUI
                              CAREA_UNIT=0, # number of cells (less accurate but faster)
                              METHOD=0     # deterministic 8
                            ))
  return(out_grid)
  }

#' @param point either SpatialPoints* or character path to points used to intersect grid. Must have same CRS as the grid
#' @param grid grid whose values will be sampled
#' @param saga.env
#' @return vector of grid values the same length as the number of points
#' @export
#' @keywords internal
SampleRasterRS <- function(point, grid, saga.env, verbose=F){
  print("Sampling raster...")
  if (grepl("Spatial", class(point))){
    rgdal::writeOGR(point, dsn=file.path(saga.env$workspace, "station.shp"), layer="station.shp", driver="ESRI Shapefile")
    point <- file.path(saga.env$workspace, "station.shp")
  }
  sampled <- gsub("\\.shp", "_gridsample\\.shp", point)
  rsaga.geoprocessor(lib='shapes_grid', module=0, env=saga.env, show.output.on.console = verbose,
                     param=list(
                       SHAPES=point,
                       GRIDS=grid,
                       RESULT=sampled,
                       INTERPOL=0
                     ))
  x <- invisible(rgdal::readOGR(sampled))
  x <- x@data[,ncol(x@data)]
  return(x)
}
