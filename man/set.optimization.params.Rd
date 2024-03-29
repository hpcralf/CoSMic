% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_functions.R
\name{set.optimization.params}
\alias{set.optimization.params}
\title{Setup of optimization parameters}
\usage{
set.optimization.params(
  opt.target.icu = TRUE,
  opt.target.deaths = FALSE,
  opt.target.region = "state",
  opt.names = c("SH07", "SH08", "SH09"),
  opt.lb = c(0, 0, 0),
  opt.ub = c(1, 1, 1),
  opt.pop.size = 20,
  opt.max.iter = 10,
  use.sug.sol = FALSE,
  opt.filter = NULL,
  ep,
  sp,
  pspc
)
}
\description{
Setup of optimization parameters
}
