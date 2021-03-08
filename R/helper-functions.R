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
#' Helps with delaying and smoothing changes in R0
#' 
#' The function smoothes numeric vectors either by logistic or linear
#' interpolation.
#'
#' @export
attenuate <- function(x,steps=5,type="logistic") {
    
    ## Find changes in x -----------------------------------------
    changex <- which(diff(x)!=0)
    nchange <- length(changex)
    
    ## Apply logistic ------------------------------------------------------
    if(type=="logistic" & nchange>0) {
        for(change in changex) {
            
            ## Get start and final value -------------------------
            start <- x[change]
            
            if ( (change+steps+1) <= length(x) ) {
                s <- steps
            } else {
                s <- length(x)-change-1
            }
            
            final <- x[change+s+1]
            
            ## Interpolation values ------------------------------
            logval <- seq(-5,5,length.out=s)
            logval <- 1 - 1/(1+exp(-logval))
            
            ## Interpolate ---------------------------------------
            newvalues <- start*logval + final *(1-logval)
            x[(change+1):(change+s)] <- newvalues
        }
    }
    
    ## Apply linear --------------------------------------------------------
    if(type=="linear" & nchange>0) {
        for(change in changex[-nchange]) {
            
            ## Get start and final value -------------------------
            start <- x[change]
            
            if ( (change+steps+1) <= length(x) ) {
                s <- steps
            } else {
                s <- length(x)-change-1
            }
            
            final <- x[change+s+1]
            
            ## Interpolation values ------------------------------
            linval <- seq(1,0,length.out=s+2)
            linval <- linval[-c(1,length(linval))]
            
            ## Interpolate ---------------------------------------
            newvalues <- start*linval + final *(1-linval)
            x[(change+1):(change+s)] <- newvalues
        }
    }
    
    return(x)
}

################################################################################
#' Extraploate 0effects beyond determined values
#' 
#' The function extrapolates R0effects based on different methods.
#'
#' @param R0effect A data.frame with R0effects per week and region.
#' @param sp An object with static CoSMic model parameters.
#' 
#' @param method Method by which to extrapolate. Supported values are:
#'               \code{"constant-weekly"}: Extrapolates constantly the R0effect
#'               of week \code{base} to the next \code{length} weeks.
#'               \code{"constant-daily"}: Determines the averaged daily R0effect
#'               from the last \code{length.days} and extraploates it constantly
#'               to the next \code{length} weeks.
#'               *Defaults to:* \code{"constant-weekly"}.
#' @param base The week based on which to extrapolate. If not given the last
#'             week i.e. dim(R0effect)\[1\] is used.
#'             *Defaults to:* \code{NULL}.
#' @param length Number of week to extrapolate after base.
#'               *Defaults to:* \code{8}.
#' @param length.days Number of days to take into account when extraploation
#'                    based on daily quantities is done.
#'                    *Defaults to:* \code{14}.
#' 
#' @export
setup.projection <- function(R0effect, sp,
                             method="constant-daily",
                             base=NULL,length=8,length.days=14) {
    
    if (is.null(base)) {
        base <- dim(R0effect)[1]
    }

     ## Extrapolate weekly value constantly --------------------------
    if (method == "constant-weekly") {
        R0effect[(base+1):(base+length),] <- R0effect[base,]
    }

    ## Constant extrapolation of averaged daily value from last length.days ----
    if (method == "constant-daily") {

        ## To allow for correct smoothing by attenuate ---------------
        R0effect[base+1,]<-R0effect[base,]
        R0effect[base+2,]<-0
      
        ## Expand weekly R0effects to days by R0change ---------------
        R0effect.d <- do.call(rbind,
                              lapply(c(1:length(sp$R0change)),
                                     function(x){R0effect[rep(x,sp$R0change[[x]][2]-
                                                                sp$R0change[[x]][1]+1),]
                                     })
                              )

        R0effect.ds <- apply(R0effect.d,2,attenuate,steps=sp$R0delay_days,type=sp$R0delay_type)

        R0.daily <- apply(R0effect.d[sp$R0change[[base]][2]:(sp$R0change[[base]][2]-14),],2,mean)
        
        R0effect[(base+1):(base+length),] <- rep(R0.daily,each=length)
    }

    return(R0effect)
    
}

################################################################################
#' Export csv
#' 
#' The function exports the contents of rr as csv files.
#'
#' @export
export.csv <- function(rr) {

    
    
}
