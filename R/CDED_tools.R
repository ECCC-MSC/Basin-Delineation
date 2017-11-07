# USGS national grid
# 1 arc-sec https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/1/ArcGrid/ (covers canada too!)


#' Get URL of CDED NTS tile
#' @param NTS character NTS string in format DDDLDD or DDDL (d=digit, L=letter) e.g 075M05 or 014K
#' @return character path to CDED DEM tile for selected NTS sheet
#' @export
#' @keywords internal
CDEDGetTilePath <- function(NTS){
  NTS <- tolower(NTS)
  resolution <- switch(as.character(nchar(NTS)), '6'='50k_dem', '4'='250k_dem')
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec"  ###  these files don't tile properly at 250k
  FTP.path <- paste(basepath, resolution, substr(NTS, 1,3), NTS, sep='/')
  FTP.path <- paste(FTP.path, '.zip', sep='')
  return(FTP.path)
}

#' Get URL of CDEM NTS tile
#' @param NTS character NTS string in format DDDLDD or DDDL (d=digit, L=letter) e.g 075M05 or 014K
#' @return character path to CDED DEM tile for selected NTS sheet
#' @export
#' @keywords internal
CDEMGetTilePath <- function(NTS){
  NTS <- toupper(NTS)
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation/cdem_mnec"
  DEM.name <- paste("cdem_dem",NTS, "tif.zip", sep="_")
  FTP.path <- paste(basepath, substr(NTS, 1,3), DEM.name, sep='/')
  return(FTP.path)
}

#' Download NTS Sheet
#' @param NTS character string, name of NTS tile (e.g. '075C' or '014D03'). If specified as
#' a 250k tile (4 characters), a 250k tile will be downloaded. If specified as a 50k tile (6 character)
#' then a 50k tile will be downloaded.
#' @param NTS.dir character string, path to directory where NTS tiles are stored. This function creates a nested
#' file hierarchy that matches that of the FTP site
#' @param replace logical, degaults to FALSE, whether or not files should be re-downloaded if they already exist
#' in the NTS.dir directory
#' @param product either 'CDEM' or 'CDED', which data product to get
#' @description downloads DEM tiles from "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation" and unzips
#' them to a directory
#' @export
DownloadSingleNTS <-  function(NTS, NTS.dir, replace=F, product='CDED'){
  output <- T
  FTP.path <- switch(product, 'CDED'=CDEDGetTilePath(NTS), 'CDEM'=CDEMGetTilePath(NTS))
  file.paths <- rev(rev(strsplit(FTP.path,'/')[[1]])[1:3])
  dir.create(file.path(NTS.dir, paste(file.paths[1:2], collapse='/')), showWarnings = F, recursive=T)
  destfile <- file.path(NTS.dir, paste(file.paths, collapse='/'))
  dest.dir <- gsub("\\.zip", "", destfile)
  if (!replace & dir.exists(dest.dir)){
    cat(sprintf("%s exists locally and was not downloaded\n\n", NTS))
  }else{
  output <- DownloadAndUnzip(url=FTP.path, destfile=destfile, exdir=dest.dir)
  }
  if (!is.null(output)){
    dem <- list.files(path=dest.dir, pattern="dem[ew_].*(\\<tif\\>|\\<dem\\>)$", full.names = T)
    return(dem)
  }else{
    c(NA,NA)
  }}


#' Download and Unzip NTS
#' @description a subroutine of \link{DownloadSingleNTS}
DownloadAndUnzip <- function(url, destfile, exdir){
  output <- tryCatch(
    {
      download.file(url = url, destfile = destfile, quiet = FALSE, mode =  "wb", cacheOK = TRUE)
      unzip(zipfile = destfile, exdir = exdir)
      return(1)
    },
    error = function(e){
      message(paste("NTS tile does not exist"))
      return(NULL)
    },
    warning = function(w){
      message(paste("NTS tile does not exist"))
      #return(c(NA,NA))
      return(NULL)
    })
}

#' Download Multiple NTS
#' @description a wrapper for \link{DownloadSingleNTS} that allows a vector of NTS sheets to be specified
#' @param NTS character string, name of NTS tile (e.g. '075C' or '014D03'). If specified as
#' a 250k tile (4 characters), a 250k tile will be downloaded. If specified as a 50k tile (6 character)
#' then a 50k tile will be downloaded.
#' @param NTS.dir character string, path to directory where NTS tiles are stored. This function creates a nested
#' file hierarchy that matches that of the FTP site.
#' @param force.50k logical, returns 50k tiles instead of 250k tiles even if the tiles are specified as 250k tiles.
#' e.g. using "075C" would use 075C01, 075C02, 075C03 and so on.
#' @param product either 'CDEM' or 'CDED'. Which dataset to download
#' @export
DownloadMultipleNTS <- function(NTS, NTS.dir, force.50k=F, product='CDED'){
  if (!all(grepl("^\\d{3}[[:alpha:]](\\d{2})?$", NTS))){
    stop("Bad format for one or more NTS strings")
  }
  if (force.50k & product!='CDEM'){
    NTS <- unlist(lapply(NTS, function(x) paste(x, sprintf("%02d", seq(1,16)), sep='')))
  }
  files <- t(sapply(NTS, DownloadSingleNTS, NTS.dir=NTS.dir, replace=F, product=product))
  files <- as.character(na.omit(files))
  return(files)
}

#' CompileNTS
#'
#' @description Downloads a set of NTS tiles and moasaicks them into one SAGA grid file.
#' @param NTS character string, name of NTS tile (e.g. '075C' or '014D03'). If specified as
#' a 250k tile (4 characters), a 250k tile will be downloaded. If specified as a 50k tile (6 character)
#' then a 50k tile will be downloaded.
#' @param NTS.dir character string, path to directory where NTS tiles are stored. This function creates a nested
#' file hierarchy that matches that of the FTP site.
#' @param output.dir character string specifying where output file should be saved
#' @param force.50k logical, returns 50k tiles instead of 250k tiles even if the tiles are specified as 250k tiles.
#' e.g. using "075C" would use 075C01, 075C02, 075C03 and so on.
#' @param product character, which DEM product to use.  One of ('CDED', 'CDEM')
#' @return name of output DEM
#' @examples
#' \dontrun{
#' CompileNTS("035D", "C:/temp", "C:/temp", "C:/temp/NTSmosaic.sdat", force.50k=T)
#' }
#' @export
CompileNTS <- function(NTS, NTS.dir, output.dir, force.50k=F, product='CDED'){
  files <- DownloadMultipleNTS(NTS=NTS, NTS.dir=NTS.dir, force.50k=force.50k, product=product)
  dstfile <- file.path(output.dir, "NTS_mosaic.sdat")
  gdal_mosaic(srcfiles = files, dstfile = dstfile, of = "SAGA",srcnodata = -32767 )
  output <- gdal_warp2SAGA(dstfile, srcnodata=-32768,
                 outputCRS="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  file.remove(dstfile)
  return(output)
}

#' OverlayNTS
#' @description Covers a shape with an NTS grid
#' @param geom1 An R Spatial* object to be 'buffered' with the DEM
#' @param tileindex spatialpointspolygon of NTS tile index or path to shapefile of tile index
#' @param NTS.dir directory where NTS tiles are stored
#' @param output.dir where output DEM should be saved
#' @param tilename Name of column in tileindex that gives NTS sheet number
#' @param tol (optional) by how many metres to buffer geom1
#' @param ... optional arguments passed to \link{CompileNTS}
#' @details Builds a DEM from CDED that either extends a fixed distance on either
#' side of a spatialobject (if a tol(erace) is specified) or intersects the NTS index (if
#' tileindex is specified)
#' @return name of output DEM
#' @export
OverlayNTS <- function(geom1, tileindex, NTS.dir, output.dir, tilename="NTS_SNRC", tol, ...){
  # get NTS index
  if (!missing(tol)){  # use geometry points
    tile.names <- unlist(lapply(rcanvec::nts(bbox=ExpandBBox(geom1, tol)), paste, collapse=''))
  }else if (!missing(tileindex)){               # use tile index
    tile.names <- TileIndex(geom1, tileindex, tilename)
  }

  # build and mosaic
  grid <- CompileNTS(NTS=tile.names, NTS.dir=NTS.dir, output.dir=output.dir, ...)

  # return grid name
  return(grid)
  }
