#===============================================================================
#0. Set-up - load in the required libraries and codes
options(stringsAsFactors = F)

library(HYDAT)
library(sp)
library(rgdal)
library(data.table)
library(RSQLite)
library(foreign)
library(rgeos)
library(raster)
library(spatialEco)
library(scales)
library(maptools)
library(RSAGA)
library(gdalUtils)
library(rcanvec)

# source files
source('M:/R/HydroTool/R/misc.R', echo=TRUE)
source('M:/R/HydroTool/R/DEM_tools.R', echo=TRUE)
source('M:/R/HydroTool/R/SAGA_Tools.R', echo=TRUE)
source('M:/R/HydroTool/R/HYDROBASINS_Rscript.R', echo=TRUE)
source('M:/R/HydroTool/R/spatial_stations.R', echo=TRUE)
source('M:/R/HydroTool/R/Geonames_db.R', echo=TRUE)
source('M:/R/HydroTool/R/station_names.R', echo=TRUE)
source('M:/R/HydroTool/R/basin_creation.R', echo=TRUE)

##############
#' We'll need to set some directory paths, make sure to use absolute paths
####

HYBAS_Limits_dir <- "C:/NB/wkdir/HB_basin" # where should intermediate HYBAS-derived basins be saved?
Final_output_dir <- "C:/NB/wkdir/final_basin" # where should final basins be saved?
DEM_basin_dir <- "C:/NB/wkdir/DEM_basin" # where should intermediate DEM-derived basins be saved?
HYDROSHEDS_DEM_path <- "C:/Users/brownn/Downloads/HydroSHEDS_DEM_3s"  # what is root directory of hydrosheds DEM tiles
CDED_DEM_path <-  "C:/NB/wkdir/DEM"  # what is root directory of CDED DEM tiles (these will be downloaded as necessary)
temp_dir <- "C:/NB/wkdir/tempfiles"  # directory for temporary files
saga_dir <- "C:/NB/saga222/"         # directory containing saga executable and libraries
lakes_dir= 'C:/NB/wkdir/lakes'       # where should lake polygons be saved?

##############
#' and also create a saga geoprocessing object (ignore the warning if it complains about the path)
####
envi <- rsaga.env(path=saga_dir, workspace=temp_dir)

##############
#' next, we need hydrobasins polygons data (here, north america level 12)
#############
HYBAS <- readOGR("C:/Users/brownn/Downloads/Canada_Clipped_Hydrobasins/hybas_na_lev12_v1c.shp", stringsAsFactors = F)

##############
#' next, we want station coordinates, we'll get them from a HYDAT sqlite database
#'
#############
sqlite.path <- "C:/Program Files (x86)/ECCC/ECDataExplorer/Database/Hydat.sqlite3"
con <- dbConnect(RSQLite::SQLite(), sqlite.path)
stations <- SpatialHydat(con)
dbDisconnect(con)

##############
#' Let's use 08NE123 as an example, it is in a hydrobasin with no upstream basins,
#' so it gets a code "T" (for (T)erminating-basin).  It is inside hydrobasin 7120262830
#' We'll make a dataframe with this information
#'
#############
stn <- stations[stations$station_number == '08NE123', ]
h_id = HYBAS$HYBAS_ID[which(gIntersects(stn, HYBAS, byid=T))] #7120262830
df <- data.frame(list(station_number = '08NE123', code="T", HYBAS_ID = h_id), stringsAsFactors = F)


##############
#' Now this can be used to create a hydrobasins delineation for that station
#'
#############
bs1 = HYBASBasinLimits(df = df, HYBAS=HYBAS, code=df$code, outdir=HYBAS_Limits_dir)

##############
#' Let's try again with a slightly larger basin. station 08HD032 is within hydrobasin
#' 7120262830 which has two more hydrobasin polygons upstream of it. The station
#' is on the main branch below the confluence, so it has a code of 'Pcb'
#'
#############
st2 <- stations[stations$station_number == "08HD032", ]
h_id2 = HYBAS$HYBAS_ID[which(gIntersects(st2, HYBAS, byid=T))] #7120262830
df2 <- data.frame(list(station_number = '08HD032', code="Pcb", HYBAS_ID = h_id2))
bs2 = HYBASBasinLimits(df = df2, HYBAS=HYBAS, code=df2$code, outdir=HYBAS_Limits_dir)

##############
#' Now, DEM-derived basins can be calculated for these two stations, using the
#' Hydrobasins boundaries to determine what DEM area to clip out of the mosaic.
#'
#' For very large basins, it will take a long time to calculate the entire basin,
#' instead, just the
#'
#############
DEMDrainageBasin(stn,
                 DEM_path = CDED_DEM_path,
                 DEM_source= "CDED",
                 saga_env = envi,
                 prelim_basin = readOGR(bs1),
                 outdir = DEM_basin_dir,
                 clip_size = 5000, # add 5 km to hydrobasin for clipping
                 pointbuffer = 30,
                 upstr_threshold = 8000) # use ~8000 for single-hybas delineations


DEMDrainageBasin(st2,
                 DEM_path = CDED_DEM_path,
                 DEM_source= "CDED",
                 saga_env = envi,
                 prelim_basin = readOGR(bs2),
                 outdir = DEM_basin_dir,
                 clip_size = 5000, # add 5 km to hydrobasin for clipping
                 pointbuffer = 30,
                 upstr_threshold = 8000) # use ~35000 for larger basins

##############
#' DEM-derived basins can be used as-is, or they can be clipped to the
#' limits of the hydrobasins-derived delineation to ensure there is no overlap.
#' For large basins where just the outlet is delineated with the DEM, this step
#' adds in the remainder of the basin
##############

CombineWatersheds(station          = '08NE123',
                  shedspoly_folder = HYBAS_Limits_dir,
                  dempoly_folder   = DEM_basin_dir,
                  out_folder       = Final_output_dir,
                  station_point    = stn)

CombineWatersheds(station          = '08HD032',
                  shedspoly_folder = HYBAS_Limits_dir,
                  dempoly_folder   = DEM_basin_dir,
                  out_folder       = Final_output_dir,
                  station_point    = st2)

##############
#' For stations on lakes, things are a bit more complicated because the
#' station won't necessarily be alongside the 'flow line'.  In these cases,
#' it is possible to use lake polygons to help delineate the catchment.
#'
#' For some stations, it is possible to get these lake polygons automatically
#' using the station name and the canadian geographic names database API.
#'
#' in this case, the station 02JB012 is on Lac Simard, and the classification code
#' is 'Ic'.  In this example, the hydrosheds DEM is used, which generally speeds
#' up processing time compared to CDED at the cost of accuracy.
##############
st3 = stations[stations$station_number == '02JB012',]
Lake_Station(station_point=st3,
              HYBAS = HYBAS,
              code =  "Ic",
              lakes_folder= lakes_dir,
              DEM_path = HYDROSHEDS_DEM_path,
              DEM_source = 'SHEDS',
              saga_env=envi,
              outdir_DEM = DEM_basin_dir,
              outdir_HYBAS = HYBAS_Limits_dir)

CombineWatersheds(station          = '02JB012',
                  shedspoly_folder = HYBAS_Limits_dir,
                  dempoly_folder   = DEM_basin_dir,
                  out_folder       = Final_output_dir,
                  station_point    = st3)


#===============================================================================
#1. classify basins (manually)


#1b subset stations
hb_complete <- list.files(HYBAS_Limits_dir, pattern = 'shp$', full.names = T)

#===============================================================================
#2. Create DEM watershed for small and large basins separately
df <- df[-c(1,2),]

for (name in df[!grepl('[Xdu]', df$code) & !grepl('[L]', df$flag), c("station_number")]){

  hyb <- df[df$station_number == name,]
  err <- F

  # get preliminary hydrobasins area
  prelim <- hb_complete[grepl(name, hb_complete)]
  if (length(prelim) == 0){
    print('no hybas area')
    next()
  }
  prelim <- prelim[1]
  prelim <- readOGR(prelim)
  print(prelim@bbox)
  print(prod(prelim@bbox[,2] - prelim@bbox[,1])) # debug
  # for 'smaller' basins, cover entire HB_basin with DEM and calculate upstream
  if (prod(prelim@bbox[,2] - prelim@bbox[,1]) < 2){

    # CDED area: 0.03 x 0.03 = 9e-4
    # SHEDS area: 0.09 x 0.09 = 0.0081
    upstr_threshold <- ifelse(nrow(prelim@data)==1, 5000, (sum(area(prelim)) * 1e-6) / (2 * 9e-4))
    DEM_source <- 'CDED'
    DEM_path <- CDED_DEM_path
    dL <- 10000
  }else{
    print('too big for now. skipping...')
    next() # SKIP FOR NOW
    # otherwise, for larger basins, just get outlet
    prelim <- prelim[which.max(prelim$UP_AREA), ]
    dx <- LongitudeLength(hyb$Latitude) * (prelim@bbox[1, 2] - prelim@bbox[1, 1]) * 1000 # in m
    dy <- (prelim@bbox[2, 2] - prelim@bbox[2, 1]) * 111.321 * 1000 # in m
    dL <- max(dx,dy)  + 5000
    prelim <- NULL

    if (grepl('[c]', hyb$code)){
      upstr_threshold <- 35000
      DEM_source <- 'CDED'
      DEM_path <- CDED_DEM_path
    }else if(grepl('[sT]', hyb$code)){
      upstr_threshold <- 8000
      DEM_source <- 'CDED'
      DEM_path <- CDED_DEM_path
    }else{
      err <- T
    }
  }


  if (!err){

    result <- try(
      DEMDrainageBasin(stations[stations$station_number==name,],
                       DEM_path = DEM_path,
                       DEM_source= DEM_source,
                       saga_env = envi,
                       prelim_basin = prelim,
                       outdir = DEM_basin_dir,
                       clip_size = dL,
                       pointbuffer = 30,
                       upstr_threshold=upstr_threshold)
    )
  }

  if(class(result)=="try-error" | err){
    writeLines(result,  file.path(DEM_basin_dir, paste0(name,"_failed.txt")))
  }
}


#===============================================================================
#2.5. Check DEM watersheds that might be suspect


#===============================================================================
#3. Create HYBAS limits

# for (station in df[!grepl('[du]', df$code) & !grepl('[L]', df$flag), c("station_number")]){
#
#   result <- try(HYBASBasinLimits(df = df[df$station_number==station,],
#                                    HYBAS=HYBAS,
#                                    code=df$code[df$station_number==station],
#                                    outdir=HYBAS_Limits_dir)
#                 )
#
#   if(class(result)=="try-error"){
#     writeLines(result,
#               file.path(HYBAS_Limits_dir,
#                         paste0(station,"_failed.txt")))
#   }else if (file.exists(file.path(HYBAS_Limits_dir,
#                                   paste0(station,"_failed.txt")))){
#   file.remove(file.path(HYBAS_Limits_dir,
#                         paste0(station,"_failed.txt")))
# }
# }


#===============================================================================
#4. Combine HYBAS and DEM basin polygons
#
# for (name in df[!grepl('[du]', df$code) & !grepl('[L]', df$flag), c("station_number")]){
#   # result <- try(
#
#     CombineWatersheds(name, HYBAS_Limits_dir,
#                                   DEM_basin_dir,
#                                   Final_output_dir,
#                                   station_point=stations[stations$station_number == name,])
#   )
#
# if(class(result)=="try-error"){
#   writeLines(result,  file.path(Final_output_dir,
#                           paste0(name,"_failed.txt")))
#
# }else if(file.exists(file.path(Final_output_dir,
#                                paste0(name,"_failed.txt")))){
#   file.remove(file.path(Final_output_dir,
#                         paste0(name,"_failed.txt")))
# }
# }

#===============================================================================
#5.  Make Lakes

#
# for (name in df[grepl('[L]', df$flag), c("station_number")]){
#   result <- try(
#     Lake_Station(
#       station_point=stations[stations$station_number == name,],
#       HYBAS = HYBAS,
#       lakes_folder= lakes_dir,
#       DEM_path = HYDROSHEDS_DEM_path,
#       DEM_source = 'SHEDS',
#       saga_env=envi,
#       outdir_DEM = DEM_basin_dir,
#       outdir_HYBAS = HYBAS_Limits_dir)
#   )
#   if(class(result)=="try-error"){
#     writeLines(result,  file.path(HYBAS_Limits_dir,
#                             paste0(name,"_failed.txt")))
# }
# }
#
# for (name in df[grepl('[L]', df$flag), c("station_number")]){
#   result <- try(
#     CombineWatersheds(name,
#                       shedspoly_folder=HYBAS_Limits_dir,
#                       dempoly_folder=DEM_basin_dir,
#                       out_folder=Final_output_dir,
#                       station_point=stations[stations$station_number == name,])
#   )
#   if(class(result)=="try-error"){
#     writeLines(result,  file.path(Final_output_dir,
#                             paste0(name,"_failed.txt")))
# }
# }
