#' German population structure
#'
#' The German population structure on county level (NUTS-3) stratified by age groups
#' and sex as off 31st of December 2018.
#'
#' @docType data
#'
#' @usage data(pop)
#'
#' @format An object of class \code{"data.frame"} with  five columns.
#' 1. \code{"dist_id"} - An integer representing the county'S unique identifier.
#' 2. \code{"date"} - Date of data publication.
#' 3. \code{"sex"} - Sex of the respective age group.
#' 4. \code{"age_gr"} - The age group.
#' 5. \code{"total"} - Inhabitants of the county in the respective age group and
#'    with the respective sex.
#'
#' @keywords datasets
#'
#' @references Federal Statistical Office of Germany. Kreisfreie Städte und
#' Landkreise nach Fläche, Bevölkerung und Bevölkerungsdichte, (2018),
#' \href{https://www.destatis.de/DE/Themen/Laender-Regionen/
#' Regionales}{https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales}.
#'
#' @source \href{https://www.destatis.de/DE/Themen/Laender-Regionen/
#' Regionales/Gemeindeverzeichnis/Administrativ/Archiv/Standardtabellen/
#' 04_KreiseVorjahr.html}{Statistisches Bundesamt} 
#' 
"pop"
