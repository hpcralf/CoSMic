% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_functions.R
\name{map.R0effects}
\alias{map.R0effects}
\title{Map R0effects from NUTS-1 to NUTS-2}
\usage{
map.R0effects(R0effect.nuts2, R0effect.states, rows = NULL)
}
\arguments{
\item{ep}{An execution parameter list as decribed in \code{\link[=set.exec.params]{set.exec.params()}}.}

\item{sp}{A list with static model parameters as described in \code{\link[=set.static.params]{set.static.params()}}.}
}
\description{
The function maps R0effects on NUTS-1 i.e. German state level to R0effects
on NUTS-2 level
}
