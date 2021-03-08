#' Structure of German counties
#'
#' The German county structure representing the NUTS-3 level for Germany and by
#' that the spatial simulation structure in the CoSMic default setup.
#'
#' @docType data
#'
#' @usage data(counties)
#'
#' @format An object of class \code{"data.frame"} with  four columns.
#' 1. \code{"dist_id"} - An integer representing the county's unique identifier.
#' 2. \code{"Name"} - The name of the county.
#' 3. \code{"Area"} - The counties area in km².
#' 4. \code{"Inhabitants"} - The population of the county.
#'
#' @keywords datasets
#'
#' @references Federal Statistical Office of Germany. Kreisfreie Städte und
#' Landkreise nach Fläche, Bevölkerung und Bevölkerungsdichte, (2018),
#' \href{https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales}
#' {https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales}.
#'
#' @source \href{https://www.destatis.de/DE/Themen/Laender-Regionen/
#' Regionales/Gemeindeverzeichnis/Administrativ/Archiv/Standardtabellen/
#' 04_KreiseVorjahr.html}{Statistisches Bundesamt} 
#'
"counties"
