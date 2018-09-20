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
  resolution <- switch(as.character(nchar(NTS)), "6" = "50k_dem", "4" = "250k_dem")
  basepath <- "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec"  ###  these files don't always tile properly at 250k

  ftp_path <- paste(basepath, resolution, substr(NTS, 1, 3), NTS, sep = "/")
  ftp_path <- paste0(ftp_path, ".zip")

  return(ftp_path)
}



#==============================================================================
#' @title Get URL of CDEM NTS tile
#'
#' @param NTS character, NTS string in format DDDLDD or DDDL (D=digit, L=letter)
#' e.g 075M05 or 014K
#'
#' @param product character string specifying which data product to get (one of
#' "CDEM', 'CDED', 'CDSM')
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
  DEM_name <- paste("cdem_dem", NTS, "tif.zip", sep = "_")
  ftp_path <- paste(basepath, substr(NTS, 1, 3), DEM_name, sep = "/")

  return(ftp_path)
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
  resolution <- switch(as.character(nchar(NTS)), "6" = "50k_dem", "4" = "250k_dem")

  # base URL for each data product
  basepath <- switch(product,
        "CDED" = "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec", #  these files don't tile properly at 250k
        "CDEM" = "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation/cdem_mnec",
        "CDSM" = "http://ftp.geogratis.gc.ca/pub/nrcan_rncan/elevation/cdsm_mnsc")

  # construct ftp path depending on which data product is to be used
  ftp_path <- switch(product,
            "CDED" = paste(basepath, resolution, substr(NTS, 1, 3), NTS, sep = "/"),
            "CDEM" = paste(basepath, substr(NTS, 1, 3),
                           paste("cdem_dem",
                           toupper(NTS), "tif", sep = "_"), sep = "/"),
            "CDSM" = paste(basepath,
                           substr(NTS, 1, 3),
                           paste0(toupper(NTS), "_cdsm_final"), sep = "/"))

  ftp_path <- paste0(ftp_path, ".zip")

  return(ftp_path)
}

#==============================================================================
#' @title Download DEM Sheet
#'
#' @param NTS character string, name of NTS tile (e.g. '075C' or '014D03') or
#' USGS NED tile. If specified as a 250k tile (4 characters), a 250k tile will
#'  be downloaded. If specified as a 50k tile (6 character) then a 50k tile will
#'  be downloaded.
#'
#' @param DEM_dir character string, path to directory where DEM tiles are stored.
#' This function creates a nested file hierarchy that matches that of the FTP site
#'
#' @param replace logical, degaults to FALSE, whether or not files should be
#' re-downloaded if they already exist in the DEM_dir directory
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
DownloadSingleDEM <-  function(DEM_id, DEM_dir, replace=F, product="CDED"){
  output <- T

  # construct ftp path to be queried
  ftp_path <- switch(product,
                     "CDED"  = GetNTSDEMTilePath(DEM_id, product = "CDED"),
                     "CDEM"  = GetNTSDEMTilePath(DEM_id, product = "CDEM"),
                     "CDSM"  = GetNTSDEMTilePath(DEM_id, product = "CDSM"),
                     "NED"   = USGSTileURL(name = DEM_id),
                     "ADEM"  = ARCTICTileURL(name = DEM_id))

  # create appropriate file name / directory structure based on ftp path
  file_paths <- rev(rev(strsplit(ftp_path, "/")[[1]])[1:3])

  dir.create(file.path(DEM_dir, paste(file_paths[1:2], collapse = "/")),
             showWarnings = F, recursive = T)

  destfile <- file.path(DEM_dir, paste(file_paths, collapse = "/"))

  if (grepl("tar\\.gz", destfile)){
    dest_dir <- gsub("\\.tar\\.gz", "", destfile)
  }else{
    dest_dir <- gsub("\\.zip", "", destfile)
  }

  # Check to see if file already exists
  if (!replace & dir.exists(dest_dir)){
    cat(sprintf("%s exists locally and was not downloaded\n\n", DEM_id))
  }else{
    output <- DownloadAndUnzip(url = ftp_path, destfile = destfile, exdir = dest_dir)
  }

  # If an appropriate file was downloaded, return the corresponding file paths
  if (output){
    pattern <- switch(product,
                      "CDED" = "dem[ew_].*(\\<tif\\>|\\<dem\\>)$",
                      "CDEM" = "dem[ew_].*(\\<tif\\>|\\<dem\\>)$",
                      "CDSM" = ".*_cdsm_final_[ew]\\.tif",
                      "NED" = "flt$",
                      "ADEM" = "reg_dem\\.tif$")

    dem <- list.files(path = dest_dir, pattern = pattern, full.names = T)

    return(dem)

  }else{
    return(c(NA, NA))
  }
  }


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
  output <- tryCatch({
    download.file(url = url, destfile = destfile,
                  quiet = FALSE, mode =  "wb", cacheOK = TRUE)

      if (grepl("tar\\.gz", destfile)){
        untar(tarfile = destfile, exdir = exdir)
      }else{
        unzip(zipfile = destfile, exdir = exdir)
      }


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
#' @param DEM_dir character string, path to directory where DEM tiles are stored.
#' This function creates a nested
#' file hierarchy that matches that of the FTP site.
#'
#' @param force50k logical, returns 50k tiles instead of 250k tiles even if the
#' tiles are specified as 250k tiles.
#' e.g. using "075C" would use 075C01, 075C02, 075C03 and so on.
#'
#' @param product either 'CDEM' or 'CDED'. Which dataset to download
#'
#' @export
#==============================================================================
DownloadMultipleDEM <- function(DEM, DEM_dir, force50k=F, product="CDED"){

  # sanity check - make sure NTS names are well-formed
  if (toupper(product) %in% c("CDED", "CDEM", "CDSM") &
        !all(grepl("^\\d{3}[[:alpha:]](\\d{2})?$", DEM))){
      stop("Bad format for one or more NTS strings")
    }

    if (force50k & product %in% c("CDED", "CDSM")){
      DEM <- unlist(lapply(DEM, function(x)
        paste0(x, sprintf("%02d", seq(1, 16)))))
    }

  # download each DEM file using the apply function
    files <- t(sapply(DEM, DownloadSingleDEM, DEM_dir = DEM_dir,
                      replace = F, product = product))

  files <- as.character(na.omit(files))
  return(files)
}

#==============================================================================
#' @title Create DEM mosaic
#'
#' @description Downloads a set of DEM tiles and moasaicks them into one SAGA
#' grid file.
#'
#' @param DEM character vector, names of NTS tiles (e.g. '075C' or '014D03') or
#'  USGS NED 1arc-sec tiles (e.g. 'n52e104'). If specified as a 250k NTS tile
#'  (4 characters), a 250k NTS tile will be downloaded. If  specified as a 50k
#'  NTS tile (6 character) then a 50k NTS tile will be downloaded.
#'
#' @param DEM_dir character string, path to directory where DEM tiles are stored.
#' This function creates a nested file hierarchy that matches that of the
#' FTP site.
#'
#' @param output_dir character string specifying where output file should be saved
#'
#' @param force50k logical, returns 50k tiles instead of 250k tiles even if the
#' tiles are specified as 250k tiles. e.g. using "075C" would use 075C01, 075C02,
#'  075C03 and so on. Ignored if using 'NED' as data product.
#'
#' @param product character, which DEM product to use.  One of ('CDED', 'CDEM',
#'  'NED')
#'
#' @return name of output DEM
#'
#' @examples
#'
#' CompileDEM("035D", "~", "~", "~/NTSmosaic.sdat", force50k = T)
#'
#' @export
#==============================================================================
CompileDEM <- function(DEM, DEM_dir, output_dir, force50k = F, product = "CDED"){
  files <- DownloadMultipleDEM(DEM = DEM, DEM_dir = DEM_dir, force50k = force50k,
                               product = product)

  dstfile <- file.path(output_dir, "NTS_mosaic.sdat")
  gdal_mosaic(srcfiles = files[file.exists(files)], dstfile = dstfile,
              of = "SAGA", srcnodata = -32767)

  output <- gdal_warp2SAGA(dstfile, srcnodata = -32768,
                 output_crs = "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96
                 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

  file.remove(dstfile)

  return(output)
}

#==============================================================================
#' @title Create DEM mosaic to cover spatial object
#'
#' @description Covers a shape with a DEM grid
#'
#' @param geom1 An R Spatial* object to be covered with the DEM
#'
#' @param tileindex spatialpointspolygon of DEM tile indices or path to
#'  shapefile of tile index
#'
#' @param DEM_dir directory where DEM tiles are stored
#'
#' @param output_dir where output DEM should be saved
#'
#' @param tileid_field Name of column in tileindex that gives DEM sheet number
#'
#' @param tol (optional) how wide should geom1 be buffered (metres) before
#' determining how big the DEM needs to be to cover it
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
OverlayDEM <- function(geom1, tileindex, DEM_dir, output_dir, tol,
                       product="CDED", tileid_field="NTS_SNRC", ...){

  # get DEM names
  if (!missing(tol)){
    # use geometry points
    if (product %in% c("CDED", "CDEM", "CDSM")){
      atscale <- ifelse(toupper(product) == "CDEM", 1, 2)
      tiles <- rcanvec::nts(bbox = ExpandBBox(geom1, tol), atscale = atscale)

      # convert to a character vector
      if (class(tiles) == "list"){
        tile_names <- sapply(tiles, paste, collapse = "")
      }else{
        tile_names <- paste(tiles, collapse = "")
      }

    }else if (product == "NED"){
      tile_names <- NEDcoverage(geom1, tol)
    }
  }else if (!missing(tileindex)){
    # use tile index
    tile_names <- TileIndex(geom1, tileindex, tileid_field)
  }

  # build and mosaic
  grid <- CompileDEM(DEM = tile_names, DEM_dir = DEM_dir, output_dir = output_dir,
                     product = product, ...)

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
USGSTileURL <- function(lon, lat, name, test=TRUE){

  baseURL <- "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1/GridFloat/"

  if (missing(name)){
    name <- USGSTileName(lon, lat, fext = ".zip")
  }else{
    name <- paste0(name, ".zip")
  }

  fullpath <- sprintf("%s%s", baseURL, name)

  if (test & http_error(fullpath)){

    if (!missing(lon)){
      name <- USGSTileName(lon, lat, fext = ".zip", v='2017')
      fullpath <- sprintf("%s%s", baseURL, name)

      }else{
      fullpath <- sub(name, paste0("USGS_NED_1_",name), fullpath)
      fullpath <- sub("\\.zip", "_GridFloat.zip", fullpath)
    }
  }

  return(fullpath)
}

#==============================================================================
#' @title USGS Tile name from coordinates
#'
#' @description Returns the base file name of a USGS NED elevation tile that
#' covers a set of coordinates
#'
#' @param lon longitude in degrees as an integer
#'
#' @param lat latitude in degrees as an integer
#'
#' @param fext file extension to append to file name
#'
#' @return a character string of format "n00w000%s"
#==============================================================================
USGSTileName <- function(lon, lat, fext="", v="2013"){
  if (v=="2013"){
    name <- sprintf("n%0.2dw%0.3d%s", abs(lat), abs(lon), fext)
  }else if (v=="2017"){
    name <- sprintf("usgs_ned_1_n%0.2dw%0.3d_gridfloat%s", abs(lat), abs(lon), fext) # July 2018 update
  }

  return(name)
}



#==============================================================================
#' @title ARCTIC DEM Tile name from coordinates
#'
#' @description Returns the base file name of an ARCTIC DEM elevation tile
#'
#' @param lon longitude in degrees
#'
#' @param lat latitude in degrees
#'
#' @param fext file extension to append to file name
#'
#' @return a character string
#==============================================================================
ARCTICTileName <- function(lon, lat, fext=""){

  # convert to NSIDC polar coords (epsg: 3413 )
  pt <- SpatialPoints(data.frame(x=lon, y=lat), proj4string = CRS("+init=epsg:4326"))
  pt <- sp::spTransform(pt, CRS("+init=epsg:3413"))
  x <- pt@coords[1, 1]
  y <- pt@coords[1, 2]

  # find which band
  ymin_band <- floor(y / 100000) * 100000
  band <- 41 + ymin_band / 100000

  # find which column
  xmin_col <- floor(x / 100000) * 100000
  col <- 41 + (xmin_col / 100000)

  # find which sub-square
  xc <- 1 + floor(x / 50000) %% 2
  yc <- 1 + floor(y / 50000) %% 2

  tile <- sprintf("%02d_%02d_%01d_%01d_5m_v2.0%s", band, col, xc, yc, fext)

  return(tile)
  }

#==============================================================================
#' @title  ARCTIC DEM Tile URL
#'
#' @description Returns the URL of a ARCTIC DEM elevation tile
#'
#' @param lon longitude in degrees as an integer
#'
#' @param lat latitude in degrees as an integer
#'
#' @param name (optional) tilename
#'
#' @return file path
#==============================================================================
ARCTICTileURL <- function(lon, lat, name){

  if (missing(name)){
    name <- ARCTICTileName(lon, lat, fext = ".tar.gz")
  }else{
    name <- paste0(name, ".tar.gz")
  }

  return(sprintf("ftp://ftp.data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/%s/%s", substr(name, 1, 5), name))
}
