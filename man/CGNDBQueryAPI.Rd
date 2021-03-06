% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geonames_db.R
\name{CGNDBQueryAPI}
\alias{CGNDBQueryAPI}
\title{Construct query for CGNDB API}
\usage{
CGNDBQueryAPI(terms, format = NULL)
}
\arguments{
\item{terms}{a list of named search terms, see details for more info}

\item{format}{Character string; one of c('json', 'html', 'csv', 'kml', 'gml',
'wkt'). Determines the format of the results returned by the API.}
}
\value{
Character string; a URL-encoded search string
}
\description{
Returns a URL that can be used to query the
Canadian Geographic Names Database using the online API
}
\details{
The terms parameter must be a list with the names corresponding to
search parameters of the CGNDB.
More details about search terms can be found at
http://www.nrcan.gc.ca/earth-sciences/geography/place-names/tools-applications/9249
Common codes include q, bbox, lat, lon, radius, category, theme, and province

Theme codes allow for only a certain category of feature to be returned. The
theme codes are:
979  Water Feature
980 	Volcanic Feature
981 	Undersea and Maritime Feature
982 	Underground Feature
983 	Unclassified
984 	Terrain Feature
985 	Populated Place
986 	Ice and Snow Feature
987 	Feature Associated with Vegetation
988 	Constructed Feature
989 	Administrative Area
}
\examples{
CGNDBQueryAPI()
CGNDBQueryAPI(list(q='Ottawa', theme=985), format='csv')
}
