#===============================================================================
#' @title Set upstream station
#'
#'
#' @description For a target station, sets a list of stations as 'upstream'
#'
#' @param station name of a target station for which upstream stations are to
#' be defined.
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @param upstream a character vector of stations in table that are upstream of
#' the target station.
#'
#' @return Does not return a value, the table is modified in-place.
#'
#' @examples
#' r <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3'))
#' set_upstream('stn1', r, 'stn2')
#'
#'  @export
#===============================================================================
set_upstream <- function(station, table, upstream, set_to=1){
  if (station %in% upstream){
    stop("A station cannot be upstream of itself")
  }
  if (length(upstream)==0){
    warning(sprintf("Station %s has no upstream stations. Matrix is unchanged",
                    station))
    return()
  }
    test <- names(table) %in% upstream
    table[station, names(table)[test] := set_to]
}

#===============================================================================
#' @title Set downstream station
#'
#'
#' @description For a target station, sets a list of stations as 'downstream'
#'
#' @param station name of a target station for which downstream stations are to
#' be defined.
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @param upstream a character vector of stations in table that are downstream of
#' the target station.
#'
#' @return Does not return a value, the table is modified in-place.
#'
#' @examples
#' r <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3'))
#' set_downstream('stn1', r, 'stn2')
#'
#'  @export
#===============================================================================
set_downstream <- function(station, table, downstream, set_to=1){
  x <- lapply(downstream, set_upstream, table=table, upstream=station, set_to=set_to)
  table #if this is missing the next call of the table will return nothing
}

#===============================================================================
#' @title Recursively set upstream stations
#'
#' @description recursively sets all stations upstream of a target station
#'
#' @param station target station to set
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @param next_upstream the next station upstream of the target station
#'
#' @examples
#' M  <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3', 'stn4'))
#' set_upstream('stn1', M, c('stn2', 'stn3'))
#' set_upstream_recursive('stn4', M, 'stn1')
#===============================================================================
set_upstream_recursive <- function(station, table, next_upstream, set_to=1){
  # get all upstream of next_upstream
  all_up <- get_upstream(next_upstream, table)

  # Set all stations upstream of next_upstream as upstream of target station
  set_upstream(station, table, c(next_upstream, all_up), set_to=set_to)
}

#===============================================================================
#' @title Recursively set downstream stations
#'
#' @description recursively sets all stations downstream of a target station
#'
#' @param station target station to set
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @param next_downstream the next station downstream of the target station
#'
#' @examples
#' M  <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3', 'stn4', 'stn5'))
#' set_upstream('stn1', M, c('stn2', 'stn3', 'stn4'))
#' set_upstream('stn2', M, 'stn3')
#' set_downstream_recursive('stn5', M, 'stn3')
#' get_downstream('stn5', M)
#===============================================================================
set_downstream_recursive <- function(station, table, next_downstream, set_to){
  # get all downstream of next_downstream
  all_down <- get_downstream(next_downstream, table)

  # Set target station as upstream of
  set_downstream(station, table = table,
                 downstream=c(all_down, next_downstream), set_to=set_to)
}

#===============================================================================
#' @title Set relative position of station
#'
#' @description, given the next-upstream and next-downstream stations, set
#' all up- and downstream stations as appropriate.
#'
#' @param station target station for which upstream-downstream relations
#' are to be set
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#' @param next_upstream (optional) the next station upstream of the target
#' station. If omitted, no upstream stations will be set.
#'
#' @param next_downstream (optional) the next station downstream of the target
#'  station. If omitted, no upstream stations will be set.
#'
#' @examples
#' M  <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3', 'stn4','stn5'))
#' set_upstream('stn1', M, c('stn2', 'stn3', 'stn4'))
#' set_upstream('stn2', M, c('stn3', 'stn4'))
#' set_upstream('stn3', M, 'stn4')
#' set_relative_position('stn5', M, next_upstream='stn3', next_downstream='stn2')
#' get_upstream('stn5', M)
#' get_downstream('stn5', M)

#===============================================================================
set_relative_position <- function(station, table,
                                     next_upstream, next_downstream){
  if (!missing(next_upstream) & !is.na(next_upstream)){
    set_upstream_recursive(station=station, table=table,
                           next_upstream = next_upstream)
  }
  if (!missing(next_downstream) & !is.na(next_downstream)){
    set_downstream_recursive(station = station, table = table,
                             next_downstream = next_downstream)
  }
}

#===============================================================================
#' @title Get upstream stations
#'
#'
#' @description return a list of stations that are upstream of a target station
#'
#' @param station target station for which a list of upstream stations is desired
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @return Does not return a value, the table is modified in-place.
#'
#' @examples
#' r <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3'))
#' set_upstream('stn1', r, 'stn2')
#' get_upstream('stn1', r)
#'
#' @export
#===============================================================================
get_upstream <- function(station, table){
 test <- which(table[station, ] == 1)
 return(names(table)[test])
}

#===============================================================================
#' @title Get downstream stations
#'
#'
#' @description return a list of stations that are downstream of a target station
#'
#' @param station target station for which a list of downstream
#' stations is desired
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @return character vector of stations upstream of target station
#'
#' @examples
#' r <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3'))
#' set_upstream('stn1', r, 'stn2')
#' get_upstream('stn1', r)
#'
#'  @export
#===============================================================================
get_downstream <- function(station, table){
  downstream <- as.logical(table[, station, with=F] == 1)
  result <- table[, 'ID', with=F][downstream]
  return(result[[1]])
}

#===============================================================================
#' @title Get a list of the most downstream stations
#'
#'
#' @description return a list of stations that are the most downstream in their
#' catchment area
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @param diagnl
#'
#' @return character vector of stations upstream of target station
#'
#' @examples
#' r <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3'))
#' set_upstream('stn1', r, 'stn2')
#' get_upstream('stn1', r)
#'
#'  @export
#===============================================================================
most_downstream <- function(table, diagnl=9){
  n_downstream <- apply(table[,-1], 2, sum)
  most_downstream <- c(F, n_downstream==diagnl)
  return(names(table)[most_downstream])
}

#===============================================================================
#' @title Create upstream matrix from station names
#'
#'
#' @description return a list of stations that are downstream of a target station
#'
#' @param stations character vector of station names
#'
#' @param diagnl integer used along the diagonal of the matrix to indicate
#' station identity (neither upstream nor downstream)
#'
#' @return an [n × n+1] data.table
#'
#' @examples
#' r <- initialize_matrix_from_stations(c('stn1', 'stn2', 'stn3'))
#'
#'  @export
#===============================================================================
initialize_matrix_from_stations <- function(stations, diagnl=9){
  M <- cbind(rep("", length(stations)), data.table(diagnl*diag(length(stations))))
  setnames(M, c("ID",stations))
  M[,"ID"] <- stations
  setkey(M, ID)

    setcolorder(M, c("ID", M[,ID]))

  return(M)
}


#===============================================================================
#' @title Add station to upstream table
#'
#' @description Creates a new table w.  Because of the limitations of the
#' data.table object on which the upstream-downstream tables are based, the object
#' is not modified in-place. Instead, a new object must be created.  A practical
#' workaround is to assign the new table object the same name as the old table
#'
#' @param stations a list of station names to be added to the table
#'
#' @param table a data.table describing upstream / downstream relationships between
#' stations.
#'
#' @examples
#' a <- initialize_matrix_from_stations(c('st1', 'st2','st8'))
#' set_upstream('st1', a, c('st2', 'st8'))
#' a <- add_station_to_table(c('st4', 'st5'), a)
#'
#' @export
#===============================================================================
add_station_to_table <- function(stations, table){
  table <- data.table::copy(table)

  #create new rows
  newstn <- initialize_matrix_from_stations(stations)
  emptyrows <- table[1:length(stations), -1]
  emptyrows[, 1:ncol(emptyrows) := 0]
  newrows <- cbind(newstn[, 1], emptyrows, newstn[, -1])

  # create new columns in existing table
  table[,paste0(stations) := 0]

  # combine rows together
  out <- rbindlist(list(table, newrows), use.names = T)
  setkey(out, "ID")
  setcolorder(out, c("ID", out[,ID]))

  return(out)
}


#===============================================================================
#'@title Check transitivity of station matrix
#'
#' @description  Finds all monitoring stations that violate the assumption
#'of transitivity for upstream/downstream relations.
#'
#' @details The transitivity assumption states that if x is upstream of y and
#'y is upstream of z, then x is upstream of z. This tests that assumption for a
#'matrix
#'
#' @param table: a data table created from initialize_matrix_from_stations. An
#'(n x n+1) data.table with first column named "ID" and all other columns
#'named according to station number. Row i and column i+1 should correspond to
#'the same station.
#'
#' @return A list with an element for every x,y station pair that violates
#' the transitivity assumption
#'
#' @export
#===============================================================================
check_transitivity <- function(table){
  M <- as.matrix(table[,-1]) # get numeric part of table
  M <- -M * (diag(nrow(M)) - 1) # set diagonal to zero

  # Find which entries in M violate transitivity assumption
  intransitive <- M != 0 | (M == 0 & ((M - M %*% M) == 0))

  #return names that correspond to FALSE values in intransitive
  out <- apply(intransitive, 1, function(x) table[["ID"]][!x])
  names(out) <- table[["ID"]]
  out <- out[lapply(out, length) > 0]
  return(out)
}

#===============================================================================
#'@title Check symmetry of station matrix
#'
#' @description  Finds all monitoring stations that violate the assumption
#'of antisymmetry for upstream/downstream relations.
#'
#' @details The antisymmetry assumption states that if x is upstream of y then
#' y is not upstream of x
#'
#' @param table: a data table created from initialize_matrix_from_stations. An
#'(n x n+1) data.table with first column named "ID" and all other columns
#'named according to station number. Row i and column i+1 should correspond to
#'the same station.
#'
#' @return A list with an element for every x,y station pair that violates
#' the antisymmetry assumption
#'
#' @export
#===============================================================================
check_symmetry <- function(table){
  M <- as.matrix(table[,-1]) # get numeric part of table
  M <- -M*(diag(nrow(M))-1) # set diagonal to zero

  # Find all entries in M where a station is both upstream and downstream of another station
  symmetric <-  (M == 1 & t(M) == 1)
  out <- apply(symmetric, 1, function(x) table[["ID"]][x])
  names(out) <- table[["ID"]]
  out <- out[lapply(out,length) > 0]
  return(out)
}

#===============================================================================
#' @title Find Upstream Stations
#'
#' @description  Finds all monitoring stations within a target basin
#'
#' @param basin
#' @param stations
#' @param target_station
#' @param stationID
#' @return A character vector of station IDs
#' @keywords internal
#' @export
#===============================================================================
upstream_stations_from_delineation <- function(basin, stations, target_station,
                                 station_number ='station_number'){
  basin <- InterpretShapefile(basin)
  basin <- SameCRS(basin, stations) # reproject basin if necessary
  #upstream.stations <- raster::intersect(stations, bind(basin,basin)) # doesn't work with just 1.. bug https://goo.gl/zUkwgJ
  upstream.stations <- sp::over(basin, stations, returnList = T)[[1]]
  upstream.stations <- upstream.stations[, station_number]
  upstream.stations <- upstream.stations[upstream.stations != target_station]
  return(upstream.stations)
}

#===============================================================================
#' @title Upstream stations from basin delineations
#'
#' @description  Use a folder of basin delineations to define upstream relations
#' and populate upstream relations matrix
#'
#' @param table a HYBAS spatial polygon object
#' @param polygon_folder The HYBAS.ID of the sub-basin from which to go upstream
#' @param station_pts spatialpointsdataframe of monitoring stations
#' @return A character vector of station IDs
#' @keywords internal
#' @export
#===============================================================================
populate_table_from_delineations <- function(table, polygon_folder,
                                                station_pts,
                                                station_number = 'station_number'){

  # get list of basin shapefiles
  basins <- list.files(polygon_folder, pattern='shp$', full.names = T)

  # iterate through table using shapefiles to find upstream stations
  for (name in names(table)[-1]){
    basin_file <- basins[grepl(name, basins)]

    if (length(basin_file) !=0){
      print(sprintf("setting upstream stations for %s", name))
      upstr <- upstream_stations_from_delineation(basin = basin_file,
                                    stations = station_pts,
                                    target_station = name,
                                    station_number = station_number)

      set_upstream(station = name,
                   table = table,
                   upstream = upstr)

    }else{
      print(sprintf("skipping station %s (cannot find basin shapefile with matching name)",
                    name))
    }
  }
}
