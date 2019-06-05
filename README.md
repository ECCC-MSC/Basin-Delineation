# README

## Introduction
These scripts were created to generate basin delineations and other value-added products from Water Survey of Canada hydrometric stations.

## Installation
download the latest build from the /built subdirectory and from your R session type:

    devtools::install_local(<path to *.tar.gz>)

## Dependencies
Most dependencies will download automatically provided that you install using `devtools::install_local` from a built `*.tar.gz` file.  Others will need to be downloaded and installed manually such as [HYDAT](https://github.com/CentreForHydrology/HYDAT)

## Code

### Components
The R package consists of the following files:

- **DEM_tools.R:** creates DEM mosaics from local files or from various online FTP servers (e.g. [CDED](http://ftp.geogratis.gc.ca/pub/nrcan_rncan/archive/elevation/geobase_cded_dnec/), [NED](https://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/1/GridFloat)) after downloading them
- **geonames_db.R:** For interacting with [canadian geographic names database](https://www.nrcan.gc.ca/earth-sciences/geography/place-names/search/9170) API
- **HYDROBASINS_Rscript.R:** For interacting with [HydroBASINS](http://www.hydrosheds.org/page/hydrobasins) tables
- **misc.R:** miscellaneous helper functions
- **SAGA_Tools.R:** Wrappers for [RSAGA](https://cran.r-project.org/web/packages/RSAGA/vignettes/RSAGA.html) functions, mostly those suited to basin delineation and DEM preprocessing
- **spatial_stations.R:** Functions to create R Spatial* objects corresponding to hydrometric stations from ECDE table exports or the [HYDAT sqlite database](https://ec.gc.ca/rhc-wsc/default.asp?n=9018B5EC-1)
- **station_names.R:** Functions to parse hydrometric station names and 

- **basin_creation.R:**  Relies on above files to create basin delineations from a few input parameters such as station names/numbers and coordinates

- **upstReam.R:** Used to construct and query large matrices describing upstream/downstream hydrometric station relationships



