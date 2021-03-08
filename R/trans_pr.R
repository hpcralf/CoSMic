#' Transition probabilities
#'
#' The dataset describes the transition probabilities to intensive care units
#' along with the chance to survive a Corona virus infction when beeing ill
#' or being ill and in intensive care stratified by age groups.
#'
#' @docType data
#'
#' @usage data(trans_pr)
#'
#' @format An object of class \code{"data.frame"} with  five columns.
#' 1. \code{"age_gr"} - Age groups with 5 year stepping from 0 - 90 years.
#' 2. \code{"sex"} - Label for stratification according t sex using labels
#'    "total", "m" and "f"
#' 3. \code{"surv_ill"} - Chance of surviving an infection.
#' 4. \code{"icu_risk"} - Risk for intensive care requirement when infected.
#' 5. \code{"surv_icu"} - Chance of surviving intensive care.
#' 
#' @keywords datasets
#'
"trans_pr"
