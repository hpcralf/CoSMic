% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_functions.R
\name{set.exec.params}
\alias{set.exec.params}
\title{Setup of execution parameters}
\usage{
set.exec.params(
  exec.procedure = "Basic-Param",
  parallel.method = "OMP",
  max.cores = 4,
  omp.cluster.dbg = FALSE,
  data.dir = "data",
  output.dir = NULL,
  model.version = "12.0",
  export_name = NULL,
  cp.write = FALSE,
  cp.time = 0,
  cp.reload = FALSE,
  cp.reload.time = 0,
  cp.dir = NULL
)
}
\arguments{
\item{exec.procedure}{Set the execution procedure. Valid values are
"Basic-Param" or "Optimization"\cr
\emph{Defaults to:} \code{"Basic-Param"}.}

\item{parallel.method}{Set the parallelization method. Valid values are
"OMP", "MPI" or "PSOCK"\cr
\emph{Defaults to:} \code{"OMP"}.}

\item{max.cores}{Set the maximum number of cores used in case
parallel.method = "OMP".\cr
\emph{Defaults to:} \code{4}}

\item{omp.cluster.dbg}{Whether std.out from workers should be captured to
a file called cl.out.\cr
\emph{Defaults to:} \code{FALSE}}

\item{data.dir}{Path to the directory from which input files are read.\cr
\emph{Defaults to:} \code{"data"}}

\item{model.version}{The model version string.\cr
\emph{Currently defaults to:} \code{12.0}}

\item{export_name}{File name addition for output files.\cr
\emph{Defaults to:} \code{<model.version>-<YYYY-MM-DD_hh:mm:ss>}}
}
\value{
A list with parameters needed to set up the execution of the CoSMic
function. The default structure is:\preformatted{$exec.procedure
[1] "Basic-Param"
$parallel.method
[1] "OMP"
$max.cores
[1] 4
$omp.cluster.dbg
[1] FALSE
$data.dir
[1] "Data"
$model.version
[1] "12.0"
$export_name
[1] "v12.0-2020-11-07_21:53:00"
}
}
\description{
Setup of execution parameters
}
