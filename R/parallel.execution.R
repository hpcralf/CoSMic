################################################################################
################################################################################
##      ___      __         _      
##     / __\___ / _\  /\/\ (_) ___ 
##    / /  / _ \\ \  /    \| |/ __|
##   / /__| (_) |\ \/ /\/\ \ | (__ 
##   \____/\___/\__/\/    \/_|\___|
##
##  COVID-19 Spatial Microsimulation  ---  For Germany  ########################
################################################################################
##
## Authors:      Christian Dudel
##               Matthias Rosenbaum-Feldbruegge
##               Sebastian Kluesener
##               Ralf Schneider
##
## Contact:      dudel@demogr.mpg.de
##               sebastian.kluesener@bib.bund.de
##
## Copyright (C) 2021 Federal Institute for Population Research (BIB),
##                    The Max Planck Institute for Demographic Research (MPIDR),
##                    High Performance Computing Center Stuttgart (HLRS)
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
################################################################################
#' Export Variables to slaves
#'
#' For convenience this function wraps the \code{exportDoMPI} and
#' \code{clusterExport} functions from the \code{doMPI} and \code{doParallel}
#' packages.
#'
#' @param ep Execution parameter list. Use [set.exec.params()] in order to
#'           create a valid layout.
#' @param cl A parallel cluster prepared by [init.parallel.execution()].
#' @param varlist Vector of character strings representing variable names to
#'                be exported to \code{cl}'s workers.
#' 
#' @export
export.to.slaves <- function (ep, cl, varlist) {

    if (ep$parallel.method == "MPI") {
        exportDoMPI(cl, varlist)
    } else {
        clusterExport(cl, varlist)
    }
}


################################################################################
#' Initalization of the parallel execution.
#'
#' The function prepares and initializes the parallel execution of the
#' [CoSMic()] model function on computer clusters in dependence from the
#' requested execution procedure and selected parallel execution method.
#'
#' @export
init.parallel.execution  <- function(ep, sp=NULL, op=NULL) {

    if ((ep$exec.procedure == "Optimization") & (is.null(op) | is.null(sp))) {
        stop("Neither sp nor op can be NULL with exec.procedure = 'Optimization'")
    }

    ## Initialize Parallel execution -------------------------------------------

    if (ep$parallel.method == "OMP") {
        ########################################################################
        ## Use shared memory parallel execution (Execution only on localhost) --
        
        library(doParallel)

        ## If std.out should be captured to file -------------------------------
        if (ep$omp.cluster.dbg) {
            cl <- makeCluster(min(detectCores(),ep$max.cores),
                          outfile="cl.out")
        } else {
            cl <- makeCluster(min(detectCores(),ep$max.cores))
        }
        
        registerDoParallel(cl)

        ## Load needed packages on slaves --------------------------------------
        ## Makes all all slaves call library(package=..., character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "dplyr",      character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "data.table", character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "rlist",      character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "lhs",        character.only = TRUE)

    }  else if (ep$parallel.method == "MPI") {
        ## ---------------------------------------------------------------------
        ## Use the message passing library for parallel communication ----------

        library(doMPI)
        cl <- startMPIcluster(verbose=TRUE)
        registerDoMPI(cl)

    }  else if (ep$parallel.method == "PSOCK") {
        ## ---------------------------------------------------------------------
        ## Use sockets for parallel communication ------------------------------

        library(parallelly)

        if (ep$exec.procedure == "Basic-Param") {
            workers <- read.csv(file="PBS_NODE.FILE",
                                stringsAsFactors=FALSE,
                                head=FALSE)[,1]
        } else {
            workers <- read.csv(file="PBS_NODE.FILE",
                                stringsAsFactors=FALSE,
                                head=FALSE)[,1]
            
            if (op$opt.pop.size > length(workers)) {
                warning(paste("Not enough hosts left in PBS-request to",
                              "calculate all indviduals in the population",
                              "in parallel. \n",
                              "Execution speed will be slowed down."))
            }
            
        }        

        ##>> Work around the ssh on Vulcan and Hawk not exporting the environment ----
        fileConn<-file("Rscript")
         
        writeLines(c("#!/bin/bash",
                     ##"export LD_LIBRARY_PATH=/opt/compiler/gnu/9.2.0/lib64",
                     ##"export LD_LIBRARY_PATH=/opt/hlrs/non-spack/compiler/gcc/9.2.0/lib:$LD_LIBRARY_PATH",
                     ##"export LD_LIBRARY_PATH=/opt/hlrs//bigdata/R/4.0.2/lib64:$LD_LIBRARY_PATH",
                     ##"export PATH=/opt/hlrs/bigdata/R/4.0.2/bin:$PATH",

                     "export LD_LIBRARY_PATH=/opt/hlrs/non-spack/compiler/gcc/9.2.0/lib64",
                     "export LD_LIBRARY_PATH=/opt/hlrs/non-spack/compiler/gcc/9.2.0/lib:$LD_LIBRARY_PATH",
                     "export LD_LIBRARY_PATH=/opt/hlrs/non-spack/bigdata/R/4.0.2/lib64:$LD_LIBRARY_PATH",
                     "export PATH=/opt/hlrs/non-spack/bigdata/R/4.0.2/bin:$PATH",

                     paste0("cd ",getwd()),
                     "Rscript $@"), fileConn)
        close(fileConn)
        system("chmod 700 Rscript")
        ##<< Work around the ssh on Vulcan and Hawk not exporting the environment ----
                
        ##if (ep$omp.cluster.dbg) {
        ##    cl <- makeCluster(workers, type = "PSOCK",
        ##                      rscript=paste(getwd(),"Rscript",sep="/"),
        ##                      outfile="cl.out")
        ##    
        ##} else {
        ##    cl <- makeCluster(workers, type = "PSOCK",
        ##                      rscript=paste(getwd(),"Rscript",sep="/"))
        ##}

        ## Since we are using R 4.0.2 backconnection of worker to master is
        ## done by slaveRSOCK() and not as from R 4.1 on by workRSOCK()
        cmd <- "parallel:::.slaveRSOCK()"

        lhost.ranks <- grep(Sys.info()["nodename"],workers)

        print(getwd())

        if (length(lhost.ranks) > 0 ) {
            cl0 <- makeClusterPSOCK(length(lhost.ranks),
                                    verbose=ep$omp.cluster.dbg,outfile = paste0(getwd(),"/cl0.out"))
        
            if ( length(workers[-lhost.ranks]) > 0 ) {
                cl1 <- makeClusterPSOCK(workers[-lhost.ranks],
                                        rscript=paste(getwd(),"Rscript",sep="/"),
                                        rscript_args=c("-e", shQuote(cmd)),
                                        verbose=ep$omp.cluster.dbg,outfile = paste0(getwd(),"/cl1.out"))
            }
        } else {
            cl1 <- makeClusterPSOCK(workers,
                                    rscript=paste(getwd(),"Rscript",sep="/"),
                                    rscript_args=c("-e", shQuote(cmd)),
                                    verbose=ep$omp.cluster.dbg,outfile = paste0(getwd(),"/cl1.out"))
        }
        
        print(length(lhost.ranks))
        print(length(workers[-lhost.ranks]))
        print(workers)

        if (exists("cl0") & exists("cl1")) {
            cl <- c(cl0,cl1)
        } else if (exists("cl0") & (!exists("cl1"))) {
            cl <- cl0
         } else if ((!exists("cl0")) & exists("cl1")) {
             cl <- cl1
         } else {
             stop("Could not start cluster")
         }
        
        
        library(doParallel)
        
        registerDoParallel(cl)
    
        ## Load needed packages on slaves --------------------------------------
        ## Makes all all slaves call library(package=..., character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "dplyr",      character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "data.table", character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "rlist",      character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "lhs",        character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "tictoc",     character.only = TRUE)
        dbg.out <- clusterCall(cl, library, package = "parallelly", character.only = TRUE)

        ## Make all slaves call setDTthreads(1) to disable ---
        ## Multithreadding in data.table processing        ---
        dbg.out <- clusterCall(cl,setDTthreads, threads=1)
        dbg.out <- clusterCall(cl,options, dplyr.summarise.inform=FALSE)

        # clusterCall(cl,load_all, "CoSMic")
        
        ## Call doParallel on workers in case of optimization --------
        ## to allow for nested parallelizm                    --------
        if (ep$exec.procedure == "Optimization") {
            dbg.out <- clusterCall(cl, library, package = "doParallel", character.only = TRUE)
        }
                
    } else {
        stop(paste("Sorry, the execution method",
                   ep$parallel.method,
                   "is not supported.",sep=" "))
    }

    ## Export main function to slaves ------------------------------------------
    export.to.slaves(exec.params, cl, "CoSMic")
    ## Export init.lhc to slaves -----------------------------------------------
    export.to.slaves(exec.params, cl, "init.lhc")
    ## Export attenuate to slaves ----------------------------------------------
    export.to.slaves(exec.params, cl, "attenuate")
    
    return(cl)
}

################################################################################
#' Finalize parallel execution execution.
#'
#' For convenience this function wraps the \code{closeCluster} and
#' \code{stopCluster} functions from the \code{doMPI} and \code{doParallel}
#' packages.
#'
#' @param ep Execution parameter list. Use [set.exec.params()] in order to
#'           create a valid layout.
#' @param cl A parallel cluster prepared by [init.parallel.execution()].
#' 
#' @export
finalize.parallel.execution  <- function(ep,cl) {

    ## Finalize Parallel execution --------------------
    if (exec.params$parallel.method == "MPI") {
        closeCluster(cl)
        mpi.finalize()    
    } else {
        stopCluster(cl)    
    }

}

