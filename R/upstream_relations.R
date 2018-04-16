library(data.table)

set_upstream <- function(station, table, upstream){
  if (station %in% upstream){
    stop("A station cannot be upstream of itself")
  }
  if (length(upstream)==0){
    warning(sprintf("Station %s has no upstream stations. Matrix is unchanged",
                    station))
    return()
  }
    test <- names(table) %in% upstream
    table[station, names(table)[test] := 1]
}


get_upstream <- function(station, table){
 test <- which(table[station, ] == 1)
 return(names(table)[test])
}


get_downstream <- function(station, table){
  downstream <- as.logical(table[, station, with=F] == 1)
  result <- table[, 'ID', with=F][downstream]
  return(result)
}

most_downstream <- function(table, diagonal_value=9){
  n_downstream <- apply(M[,-1], 2, sum)
  most_downstream <- c(F, n_downstream==diagonal_value)
  return(names(M)[most_downstream])
}


initialize_matrix_from_stations <- function(stations, diagnl=9){
  M <- cbind(rep("", length(stations)), data.table(diagnl*diag(length(stations))))
  setnames(M, c("ID",stations))
  M[,"ID"] <- stations
  setkey(M, ID)
  setcolorder(M, c("ID", M[,ID]))
  return(M)
}

add_station_to_matrix <- function(stations, table){

}
#' Check transitivity of station matrix
#'
#'@description  Finds all monitoring stations that violate the assumption
#'of transitivity for upstream/downstream relations.
#'@details The transitivity assumption states that if x is upstream of y and
#'y is upstream of z, then x is upstream of z. This tests that assumption for a
#'matrix
#'@param table: a data table created from initialize_matrix_from_stations. An
#'(n x n+1) data.table with first column named "ID" and all other columns
#'named according to station number. Row i and column i+1 should correspond to
#'the same station.
#'@return A list with an element for every
#'@export
check_transitivity <- function(table){
  M <- as.matrix(table[,-1]) # get numeric part of table
  M <- -M*(diag(nrow(M))-1) # set diagonal to zero
  intransitive <- M != 0 | (M == 0 & ((M - M%*%M) ==0))
  out <- apply(intransitive, 1, function(x) table[["ID"]][!x])
  names(out) <- table[["ID"]]
  out <- out[lapply(out,length)>0]
  return(out)
}

check_symmetry <- function(table){
  M <- as.matrix(table[,-1]) # get numeric part of table
  M <- -M*(diag(nrow(M))-1) # set diagonal to zero
  symmetric <-  (M == 1 & t(M) == 1)
  out <- apply(symmetric, 1, function(x) table[["ID"]][x])
  names(out) <- table[["ID"]]
  out <- out[lapply(out,length)>0]
  return(out)
}

#' Find Upstream Stations
#'
#'@description  Finds all monitoring stations within a target basin
#'
#'@param HYBAS: a HYBAS spatial polygon object
#'@param HYBAS.ID: The HYBAS.ID of the sub-basin from which to go upstream
#'@return A character vector of station IDs
#'@keywords internal
#'@export
FindUpstreamStations <- function(basin, stations, target_station,
                                 stationID ='station_number'){
  basin <- InterpretShapefile(basin)
  basin <- SameCRS(basin, stations) # reproject basin if necessary
  #upstream.stations <- raster::intersect(stations, bind(basin,basin)) # doesn't work with just 1.. bug https://goo.gl/zUkwgJ
  upstream.stations <- sp::over(basin, stations, returnList = T)[[1]]
  upstream.stations <- upstream.stations[,stationID]
  upstream.stations <- upstream.stations[upstream.stations != target_station]
  return(upstream.stations)
}


