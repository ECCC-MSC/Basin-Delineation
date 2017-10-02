


#' @param NTS character NTS string in format DDDLDD or DDDL (d=digit, L=letter) e.g 075M05 or 014K
#' @return character path to CDED DEM tile for selected NTS sheet
CDEDGetTilePath <- function(NTS){
  NTS <- tolower(NTS)
  resolution <- switch(as.character(nchar(NTS)), '6'='50k_dem', '4'='250k_dem')
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec"  ###  these files don't tile properly
  FTP.path <- paste(basepath, resolution, substr(NTS, 1,3), NTS, sep='/')
  FTP.path <- paste(FTP.path, '.zip', sep='')
  return(FTP.path)
}  

#' @description downloads and unzips CDED DEM tiles to a directory
CDEDDownloadTile <-  function(NTS, NTS.dir, replace=F){
  output <- T
  FTP.path <- CDEDGetTilePath(NTS)
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
    dem <- list.files(path=dest.dir, pattern="dem[we]", full.names = T)
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

#' @param force.50k logical, returns 50k tiles instead of 250k tiles even if the tiles are specified as 250k tiles.
#' e.g. using "075C" would use 075C01, 075C02, 075C03 and so on.
CompileNTS <- function(NTS, NTS.dir, output.dir, dstfile="C:/NB/NTSmosaic.sdat", force.50k=F){
  if (force.50k){
    NTS <- unlist(lapply(NTS, function(x) paste(x, sprintf("%02d", seq(1,16)), sep='')))
    }
  files <- t(sapply(NTS, CDEDDownloadTile, NTS.dir=NTS.dir, replace=F))
  files <- as.character(na.omit(files))
  file.remove(file.path(output.dir,"mosaic.sdat"),showWarnings = F)
  
  gdal_mosaic(srcfiles = files, dstfile = dstfile, of = "SAGA",srcnodata = -32767 )
  
  gdal_warp2SAGA(dstfile, srcnodata=-32768, 
                 outputCRS="+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
  return(dstfile)
}

