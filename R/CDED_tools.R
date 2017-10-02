


#' @param NTS character NTS string in format DDDLDD or DDDL (d=digit, L=letter) e.g 075M05 or 014K
#' @return character path to CDED DEM tile for selected NTS sheet
CDEDGetTilePath <- function(NTS){
  NTS <- tolower(NTS)
  resolution <- switch(as.character(nchar(NTS)), '6'='50k_dem', '4'='250k_dem')
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec"  ###  these files don't tile properly at 250k
  FTP.path <- paste(basepath, resolution, substr(NTS, 1,3), NTS, sep='/')
  FTP.path <- paste(FTP.path, '.zip', sep='')
  return(FTP.path)
}

#' @param NTS character NTS string in format DDDLDD or DDDL (d=digit, L=letter) e.g 075M05 or 014K
#' @return character path to CDED DEM tile for selected NTS sheet
CDEMGetTilePath <- function(NTS){
  NTS <- toupper(NTS)
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation/cdem_mnec"
  DEM.name <- paste("cdem_dem",NTS, "tif.zip", sep="_")
  FTP.path <- paste(basepath, substr(NTS, 1,3), DEM.name, sep='/')
  return(FTP.path)
}


#' @param NTS name of NTS tile
#' @description downloads and unzips DEM tiles to a directory
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

#' Wrapper for single download
DownloadMultipleNTS <- function(NTS, NTS.dir, force.50k=F, product){
  if (!all(grepl("^\\d{3}[[:alpha:]]{2}(\\d{2})?$"),NTS)){
    stop("Bad format for NTS string")
  }
  if (force.50k & product!='CDEM'){
    NTS <- unlist(lapply(NTS, function(x) paste(x, sprintf("%02d", seq(1,16)), sep='')))
  }
  files <- t(sapply(NTS, DownloadSingleNTS, NTS.dir=NTS.dir, replace=F, product=product))
  files <- as.character(na.omit(files))
  return(files)
}

#' @description
#' @param NTS
#' @param NTS.dir
#' @param output.dir
#' @param dstfile
#' @param force.50k logical, returns 50k tiles instead of 250k tiles even if the tiles are specified as 250k tiles.
#' e.g. using "075C" would use 075C01, 075C02, 075C03 and so on.
#' @param product character, which DEM product to use.  One of ('CDED', 'CDEM')
#' @example
#' CompileNTS("035D", "C:/temp", "C:/temp", "C:/temp/NTSmosaic.sdat", force.50k=T)
CompileNTS <- function(NTS, NTS.dir, output.dir, dstfile, force.50k=F, product='CDED'){
  files <- DownloadMultipleNTS(NTS, NTS.dir=NTS.dir, force.50k=force.50k, product=product)
  gdal_mosaic(srcfiles = files, dstfile = dstfile, of = "SAGA",srcnodata = -32767 )

  gdal_warp2SAGA(dstfile, srcnodata=-32768,
                 outputCRS="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  return(dstfile)
}

