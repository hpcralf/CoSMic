#' R0effect
#'
#' mu values for the German NUTS-2 regions representing the R0 reduction factor per
#' week and region as described by \[Klüsener-2020\]. The dataset contains mu values
#' for the simulation of all German NUTS-2 regions for 20 weeks beginning 9th of
#' March 2020.
#'
#' @docType data
#'
#' @usage data(R0effect)
#'
#' @format An object of class \code{"data.frame"} with 38 columns, one per NUTS-2
#' region, and 20 rows, one per week. If the simulation timeframe is to be extended,
#' one row per week has to be added.
#'
#' @keywords datasets
#' 
#' @references \[Klüsener-2020\] Klüsener S. et.al, Forecasting intensive care unit demand during the
#' COVID-19 pandemic: A spatial age-structured microsimulation model, (2020),
#' medRxiv,\cr
#' doi:10.1101/2020.12.23.20248761,
#' \href{https://www.medrxiv.org/content/10.1101/2020.12.23.20248761v1}{https://
#' www.medrxiv.org/content/10.1101/2020.12.23.20248761v1} .
#' 
"R0effect"
