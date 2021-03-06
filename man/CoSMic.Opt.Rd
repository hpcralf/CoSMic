% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoSMic.R
\name{CoSMic.Opt}
\alias{CoSMic.Opt}
\title{Application of the GA algorithm to CoSMic.}
\usage{
CoSMic.Opt(ep, sp, iol, pspace, sim.struc, op, cl)
}
\arguments{
\item{ep}{Execution parameter list. Use \code{\link[=set.exec.params]{set.exec.params()}} in order to
create a valid layout.}

\item{sp}{List with static model parameters. Use \code{\link[=set.static.params]{set.static.params()}} to
create a valid layout.}

\item{iol}{Input data list. Use \code{\link[=load.input]{load.input()}} to load needed fies and
\code{\link[=init.connectivity]{init.connectivity()}} in order to create a valid date layout.}

\item{pspace}{List holding the parameter space with potentially variable
model parameters. Use the setter function \code{\link[=set.pspace]{set.pspace()}} to
add parameters.}

\item{sim.struc}{List with population data. Use \code{\link[=init.spatial.population]{init.spatial.population()}}
in order to create a valid layout.}

\item{op}{List with steering parameters for the optimization process. \cr
Use \code{\link[=set.optimization.params]{set.optimization.params()}} in order to create a valid layout
and\cr \code{\link[=init.reference.data]{init.reference.data()}} in order to init the optimization
targets based on observed data.}

\item{cl}{A parallel cluster prepared by \code{\link[=init.parallel.execution]{init.parallel.execution()}}.}
}
\description{
The function applies the GA algorithm to the CoSMic simulation model
function. it uses the wrapper function \code{\link[=ff]{ff()}} as the objective function and
\code{\link[=GA.Monitor]{GA.Monitor()}} to return intermediate results during the course of the
optimization.
}
