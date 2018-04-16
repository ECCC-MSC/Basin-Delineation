#==============================================================================
#' @title Get URL of CDED NTS tile
#'
#' @param NTS character NTS string in format DDDLDD or DDDL (d=digit, L=letter)
#'  e.g 075M05 or 014K
#'
#' @return character path to CDED DEM tile for selected NTS sheet
#'
#' @export
#'
#' @keywords internal
#===============================================================================
CDEDGetTilePath <- function(NTS){
  NTS <- tolower(NTS)
  resolution <- switch(as.character(nchar(NTS)), '6'='50k_dem', '4'='250k_dem')
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec"  ###  these files don't tile properly at 250k
  FTP.path <- paste(basepath, resolution, substr(NTS, 1,3), NTS, sep='/')
  FTP.path <- paste(FTP.path, '.zip', sep='')
  return(FTP.path)
}

#==============================================================================
#' @title Get URL of CDEM NTS tile
#'
#' @param NTS character, NTS string in format DDDLDD or DDDL (d=digit, L=letter)
#' e.g 075M05 or 014K
#'
#' @param product character string specifying which data product to get (one of
#' 'CDEM', 'CDED', 'CDSM')
#'
#' @return character path to CDED DEM tile for selected NTS sheet
#'
#' @export
#'
#' @keywords internal
#==============================================================================
CDEMGetTilePath <- function(NTS, product="CDEM"){
  NTS <- toupper(NTS)
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation/cdem_mnec"
  DEM.name <- paste("cdem_dem",NTS, "tif.zip", sep="_")
  FTP.path <- paste(basepath, substr(NTS, 1,3), DEM.name, sep='/')
  return(FTP.path)
}

#==============================================================================
#' @title Get FTP path for NTS DEM tile
#'
#' @param NTS character, NTS tile name in format DDDLDD or DDDL (d=digit,
#'  L=letter) e.g 075M05 or 014K
#'
#' @param product character string specifying which data product to get (one of
#' 'CDEM', 'CDED', 'CDSM')
#'
#' @return character, FTP path to desired NTS tile
#==============================================================================
GetNTSDEMTilePath <- function(NTS, product="CDED"){
  NTS <- tolower(NTS)
  resolution <- switch(as.character(nchar(NTS)), '6'='50k_dem', '4'='250k_dem')
  basepath <- switch(product,
        "CDED" = "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec",#  these files don't tile properly at 250k
        "CDEM" = "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation/cdem_mnec",
        "CDSM" = "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation/cdsm_mnsc")
  FTP.path <- switch(product,
            "CDED" = paste(basepath, resolution, substr(NTS, 1,3), NTS, sep='/'),
            "CDEM" = paste(basepath, substr(NTS, 1,3),
                           paste("cdem_dem",
                           toupper(NTS), "tif", sep="_"), sep='/'),
            "CDSM" = paste(basepath,
                           substr(NTS, 1,3),
                           paste(toupper(NTS), "_cdsm_final", sep=""), sep='/'))
  FTP.path <- paste(FTP.path, '.zip', sep='')
  return(FTP.path)
}

#==============================================================================
#' @title Download DEM Sheet
#'
#' @param NTS character string, name of NTS tile (e.g. '075C' or '014D03') or
#' USGS NED tile. If specified as a 250k tile (4 characters), a 250k tile will
#'  be downloaded. If specified as a 50k tile (6 character) then a 50k tile will
#'  be downloaded.
#'
#' @param DEM.dir character string, path to directory where DEM tiles are stored.
#' This function creates a nested file hierarchy that matches that of the FTP site
#'
#' @param replace logical, degaults to FALSE, whether or not files should be
#' re-downloaded if they already exist in the DEM.dir directory
#'
#' @param product character string specifying which data product to get (one of
#'  'CDEM', 'CDED', 'CDSM', or 'NED')
#'
#' @description downloads DEM tiles from
#' "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation" and unzips
#' them to a directory
#'
#' @details CDED and CDEM are from NRCAN. NED uses USGS National Elevation Dataset.
#'
#' @return a vector of file paths to downloaded DEM files
#'
#' @export
#==============================================================================
DownloadSingleDEM <-  function(DEM.ID, DEM.dir, replace=F, product='CDED'){
  output <- T
  FTP.path <- switch(product,
                     'CDED'=GetNTSDEMTilePath(DEM.ID, product='CDED'),
                     'CDEM'=GetNTSDEMTilePath(DEM.ID, product='CDEM'),
                     'CDSM'=GetNTSDEMTilePath(DEM.ID, product='CDSM'),
                     'NED'=USGSTileURL(name=DEM.ID))
  file.paths <- rev(rev(strsplit(FTP.path,'/')[[1]])[1:3])
  dir.create(file.path(DEM.dir, paste(file.paths[1:2], collapse='/')),
             showWarnings = F, recursive=T)
  destfile <- file.path(DEM.dir, paste(file.paths, collapse='/'))
  dest.dir <- gsub("\\.zip", "", destfile)
  if (!replace & dir.exists(dest.dir)){
    cat(sprintf("%s exists locally and was not downloaded\n\n", DEM.ID))
  }else{
    output <- DownloadAndUnzip(url=FTP.path, destfile=destfile, exdir=dest.dir)
  }
  if (!is.null(output)){
    pattern <- switch(product,
                      'CDED'="dem[ew_].*(\\<tif\\>|\\<dem\\>)$",
                      'CDEM'="dem[ew_].*(\\<tif\\>|\\<dem\\>)$",
                      'CDSM'=".*_cdsm_final_[ew]\\.tif",
                      'NED'="^float.*flt$")
    dem <- list.files(path=dest.dir, pattern=pattern, full.names = T)
    return(dem)
  }else{
    c(NA,NA)
  }}


#==============================================================================
#' @title Download and unzip a file
#'
#' @description Downloads a zipfile from a URL and unzips it to a target
#' directory.
#'
#' @details Runs as a subroutine of \link{DownloadSingleDEM}
#'
#' @param url character url path
#'
#' @param destfile character, filepath of output zipfile
#'
#' @param exdir character,  the directory to which files are extracted
#'
#' @param rmzip logical, whether or not to remove zipfile after extraction.
#' Defaults to TRUE.
#==============================================================================
DownloadAndUnzip <- function(url, destfile, exdir, rmzip=T){
  output <- tryCatch(
    {
      download.file(url = url, destfile = destfile, quiet = FALSE, mode =  "wb",
                    cacheOK = TRUE)
      unzip(zipfile = destfile, exdir = exdir)
      if (rmzip){
        file.remove(destfile)
      }
      return(TRUE)
    },
    error = function(e){
      message(paste("Fauly DEM path or NTS tile does not exist"))
      return(FALSE)
    },
    warning = function(w){
      message(paste("Fauly DEM path or NTS tile does not exist"))
      return(FALSE)
    })
}


#==============================================================================
#' @title Download Multiple DEM files
#'
#' @description a wrapper for \link{DownloadSingleDEM} that allows a vector of
#' NTS sheets
#' or USGS DEM tile names to be specified
#'
#' @param DEM character string, name of DEM tile (e.g. '075C' or '014D03'). If
#' specified as
#' a 250k tile (4 characters), a 250k tile will be downloaded. If specified as a
#' 50k tile (6 character)
#' then a 50k tile will be downloaded.
#'
#' @param DEM.dir character string, path to directory where DEM tiles are stored.
#' This function creates a nested
#' file hierarchy that matches that of the FTP site.
#'
#' @param force.50k logical, returns 50k tiles instead of 250k tiles even if the
#' tiles are specified as 250k tiles.
#' e.g. using "075C" would use 075C01, 075C02, 075C03 and so on.
#'
#' @param product either 'CDEM' or 'CDED'. Which dataset to download
#'
#' @export
#==============================================================================
DownloadMultipleDEM <- function(DEM, DEM.dir, force.50k=F, product='CDED'){

  if (toupper(product) %in% c('CDED', 'CDEM', 'CDSM') &
        !all(grepl("^\\d{3}[[:alpha:]](\\d{2})?$", DEM))){
      stop("Bad format for one or more NTS strings")
    }

    if (force.50k & product %in% c('CDED', 'CDSM')){
      DEM <- unlist(lapply(DEM, function(x)
        paste(x, sprintf("%02d", seq(1,16)), sep='')))
    }

    files <- t(sapply(DEM, DownloadSingleDEM, DEM.dir=DEM.dir,
                      replace=F, product=product))
  files <- as.character(na.omit(files))
  return(files)
}

#==============================================================================
#' @title Create DEM mosaic
#'
#' @description Downloads a set of DEM tiles and moasaicks them into one SAGA
#' grid file.
#'
#' @param DEM character string, name of NTS tile (e.g. '075C' or '014D03') or
#'  USGS NED 1arc-sec tile (e.g. 'n52e104'). If specified as a 250k NTS tile
#'  (4 characters), a 250k NTS tile will be downloaded. If  specified as a 50k
#'  NTS tile (6 character) then a 50k NTS tile will be downloaded.
#'
#' @param DEM.dir character string, path to directory where DEM tiles are stored.
#' This function creates a nested file hierarchy that matches that of the
#' FTP site.
#'
#' @param output.dir character string specifying where output file should be saved
#'
#' @param force.50k logical, returns 50k tiles instead of 250k tiles even if the
#' tiles are specified as 250k tiles. e.g. using "075C" would use 075C01, 075C02,
#'  075C03 and so on. Ignored if using 'NED' as data product.
#'
#' @param product character, which DEM product to use.  One of ('CDED', 'CDEM',
#'  'NED')
#'
#' @return name of output DEM
#'
#' @examples
#' \dontrun{
#' CompileDEM("035D", "C:/temp", "C:/temp", "C:/temp/NTSmosaic.sdat", force.50k=T)
#' }
#'
#' @export
#==============================================================================
CompileDEM <- function(DEM, DEM.dir, output.dir, force.50k=F, product='CDED'){
  files <- DownloadMultipleDEM(DEM=DEM, DEM.dir=DEM.dir, force.50k=force.50k,
                               product=product)
  dstfile <- file.path(output.dir, "NTS_mosaic.sdat")
  gdal_mosaic(srcfiles = files, dstfile = dstfile, of = "SAGA",srcnodata = -32767)
  output <- gdal_warp2SAGA(dstfile, srcnodata=-32768,
                 outputCRS="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96
                 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  file.remove(dstfile)
  return(output)
}

#==============================================================================
#' @title Create DEM mosaic to cover spatial object
#'
#' @description Covers a shape with a DEM grid
#'
#' @param geom1 An R Spatial* object to be 'buffered' with the DEM
#'
#' @param tileindex spatialpointspolygon of DEM tile indices or path to
#'  shapefile of tile index
#'
#' @param DEM.dir directory where DEM tiles are stored
#'
#' @param output.dir where output DEM should be saved
#'
#' @param tile.id.field Name of column in tileindex that gives DEM sheet number
#'
#' @param tol (optional) by how many metres to buffer geom1
#'
#' @param ... optional arguments passed to \link{CompileDEM}
#'
#' @details Builds a DEM from elevation dataset that either extends a fixed
#' distance on either side of a spatialobject (if a tol(erace) is specified) or
#' intersects the DEM index (if tileindex is specified)
#'
#' @return name of output DEM
#'
#' @export
#==============================================================================
OverlayDEM <- function(geom1, tileindex, DEM.dir, output.dir, tol,
                       product='CDED', tile.id.field="NTS_SNRC", ...){

  # get DEM names
  if (!missing(tol)){  # use geometry points
    if (product %in% c('CDED', 'CDEM', 'CDSM')){
      atscale = ifelse(toupper(product) == 'CDEM', 1, 2)
      tiles <- rcanvec::nts(bbox=ExpandBBox(geom1, tol), atscale=atscale)
      if (class(tiles)=='list'){
        tile.names <- unlist(lapply(tiles, paste, collapse=''))
      }else{
        tile.names <- paste(tiles, collapse='')
      }
      # tile.names <- ifelse(class(tiles)=='list',
      #                      unlist(lapply(tiles, paste, collapse='')),
      #                      paste(tiles, collapse=''))
    }else if (product=='NED'){
      tile.names <- NEDcoverage(geom1, tol)
    }
  }else if (!missing(tileindex)){ # use tile index
    tile.names <- TileIndex(geom1, tileindex, tile.id.field)
  }

  # build and mosaic
  grid <- CompileDEM(DEM=tile.names, DEM.dir=DEM.dir, output.dir=output.dir,
                     product=product, ...)

  # return grid name
  return(grid)
}

#==============================================================================
#' @title USGS Tile URL
#'
#' @description Returns the URL of a USGS NED elevation tile
#'
#' @param lon longitude in degrees as an integer
#'
#' @param lat latitude in degrees as an integer
#'
#' @param name (optional)
#'
#' @return file path
#==============================================================================
USGSTileURL <- function(lon, lat, name){
  if (missing(name)){
    name <- USGSTileName(lon,lat, fext='.zip')
  }else{
    name <- paste(name, ".zip", sep='')
  }
  return(sprintf("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/GridFloat/%s", name))
}

#==============================================================================
#' @title USGS Tile name from coordinates
#'
#' @description Returns the base file name of a USGS NED elevation tile
#'
#' @param lon longitude in degrees as an integer
#'
#' @param lat latitude in degrees as an integer
#'
#' @param fext file extension to append to file name
#'
#' @return a character string of format "n00w000%s"
#==============================================================================
USGSTileName <- function(lon,lat, fext=''){
  return(sprintf("n%0.2dw%0.3d%s", abs(lat), abs(lon), fext))
}
