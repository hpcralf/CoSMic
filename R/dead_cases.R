#' Dead cases
#'
#' Dead cases for seeding during model startup.
#'
#' @docType data
#'
#' @usage data(seed_dea)
#'
#' @format An object of class \code{"data.frame"} with  ... columns.
#' 1. \code{"dist_id"} - An integer representing the county's unique identifier.
#' 2. \code{"Name"} - The name of the county.
#' 3. \code{"Area"} - The counties area in km².
#' 4. \code{"Inhabitants"} - The population of the county.
#'
#' @keywords datasets
#'
#' @references Robert Koch-Institute, COVID-19 Dashboard, (2020),\cr
#' \href{https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/nCoV_node.html}{https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/nCoV_node.html}.
#'
#' @source \href{https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/
#' dd4580c810204019a7b8eb3e0b329dd6_0}{COVID-19 Datenhub}
#'
"seed_dea"
