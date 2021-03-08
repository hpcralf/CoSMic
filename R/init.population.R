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
##==============================================================================
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
#' Initialization of the population and its spatial structure.
#'
#' The function initializes the population and its spatial structure according
#' to the layout requested by \code{sp$sim.regions}, i.e. the regions selected 
#' either at county or state level to be simulated.
#'
#' @param iol Input data list. Use [load.input()] to load needed fies and
#'            [init.connectivity()] in order to create a valid date layout.
#' @param sp List with static model parameters. Use [set.static.params()] to 
#'           create a valid layout.
#'
#' @return A list with population data. Structured as follows:
#' ```
#'  sim.struc : List of 3
#'            $ pop     : data.frame
#'                      $ dist_id: int
#'                      $ date   : chr  
#'                      $ sex    : chr
#'                      $ age_gr : int
#'                      $ total  : int
#'            $ counties: int 
#'            $ states  : data.frame
#'                      $ Code       : int
#'                      $ inhabitants: int
#'                      $ Shortcut   : chr
#'                      $ Name       : chr
#' ```
#' @export
init.spatial.population  <- function (iol,sp) {
    
    ## List of all county codes from population --------------------------------
    all_counties <- unique(iol[["pop"]]$dist_id)
    
    if( ! sp$restrict ) {
        ## ---------------------------------------------------------------------
        ## If the spacial structure is not restricted --------------------------
        
        counties <- all_counties
        states   <- iol[["states"]]
        pop      <- iol[["pop"]]

    } else {
        ## ---------------------------------------------------------------------
        ## If the spacial structure is restricted ------------------------------

        if( all(sp$sim.regions=="NONE") ) {
            ## -----------------------------------------------------------------
            ## If no spatial structure is requested ----------------------------
        
            counties <- "NONE"
            states   <- "NONE"
            
            pop         <- iol[["pop"]]
            pop$dist_id <- all_counties[1]
            pop         <- aggregate(total~dist_id+date+sex+age_gr,data=pop,sum)
        
        } else {
            ## -----------------------------------------------------------------
            ## If restricted spatial structure is requested --------------------

            if ( class(sp$sim.regions) == "character" ) {
                ## -------------------------------------------------------------
                ## if sim.regions are given by state names ---------------------
                
                ## Restrict states -----------------------------------
                states <- iol[["states"]][iol[["states"]]$Name %in%
                                          sp$sim.regions,]
                ## Restrict counties ---------------------------------
                counties <- all_counties[as.integer(all_counties/1000) %in%
                                         states$Code]

            } else {
                ## -------------------------------------------------------------
                ## if sim.regions are given by district id_id ------------------
                
                ## Restrict states -----------------------------------
                states <- iol[["states"]][iol[["states"]]$Code %in%
                                          as.integer(counties/1000),]

                ## Set counties --------------------------------------
                counties <- sp$sim.regions
            }

            
            pop <- subset(iol[["pop"]],dist_id %in% counties)
        }
        
        sim.struc <- list(pop=pop,counties=counties,states=states)

    }

    return(sim.struc)

}

################################################################################
#' Function to save the current list of execution parameters.
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' @param sim.struc List with population data. Use [init.spatial.population()]
#'                  in order to create a valid layout.
#' 
#' @export
save.spatial.population <- function(ep, sim.struc) {
    saveRDS(file=paste(ep$output.dir,
                    "spatial.population-",ep$export_name,".RDS",sep=""),
            object=sim.struc)
}

################################################################################
#' Initialize regional connectivity
#'
#' The function initializes the regional connectivity matrix according to the
#' requested regions to simulate.
#'
#' @param iol Input data list. Use [load.input()] to load needed fies and
#'            [init.connectivity()] in order to create a valid date layout.
#' @param sp List with static model parameters. Use [set.static.params()] to 
#'           create a valid layout.
#' @param ss List with population data. Use [init.spatial.population()] in 
#'           order to create a valid layout.
#'
#' @return An input data list with modified \code{connect_work} and
#'         \code{connect_total} components. See [load.input()] about details on
#'         how the input data list has to be strucutred in order to be
#'         correctly modified by this function.
#'
#' @export
init.connectivity  <- function (iol, sp, ss) {
    
    ## List of all county codes from population --------------------------------
    all_counties <- unique(iol[["pop"]]$dist_id)
    
    ## Connectivity data #######################################################

    ## Remove unneeded variable
    iol[["connect_total"]] <- subset(iol[["connect_total"]],select= -dist_id)
    
    ## Rename using only county codes
    iol[["connect_total"]] <- as.matrix(iol[["connect_total"]])
    rownames(iol[["connect_total"]]) <- paste(all_counties)
    colnames(iol[["connect_total"]]) <- paste(all_counties)

    ## Remove unneeded variable
    iol[["connect_work"]] <- subset(iol[["connect_work"]],select= -dist_id)
    iol[["connect_work"]] <- as.matrix(iol[["connect_work"]])
    
    ## Rename using county codes
    iol[["connect_work"]] <- as.matrix(iol[["connect_work"]])
    rownames(iol[["connect_work"]]) <- paste(all_counties)
    colnames(iol[["connect_work"]]) <- paste(all_counties)

    ## Restrict
    if( sp$restrict ) {
      
        ## Rows/columns to keep
        ## Keep "fake" county if no county structure
        if ( all(ss$counties=="NONE") ) {
            cols <- paste(unique(iol[["pop"]]$dist_id)[1])
        } else {
            cols <- paste(ss$counties)
        }
      
        rows <- cols
      
        ## Number of rows/cols
        ncr <- length(cols)
      
        ## Only keep rows/columns of interest
        iol[["connect_total"]] <- matrix(iol[["connect_total"]][rows,cols],nrow=ncr)
        iol[["connect_work"]]  <- matrix(iol[["connect_work"]][rows,cols],ncol=ncr)
        
        ## Rescale: Rows sum to 1
        iol[["connect_total"]] <- apply(iol[["connect_total"]],1,function(x) x/sum(x))
        iol[["connect_work"]]  <- apply(iol[["connect_work"]],1,function(x) x/sum(x))
        
        ## Transpose back to row-orientation
        iol[["connect_total"]] <- t(iol[["connect_total"]])
        iol[["connect_work"]]  <- t(iol[["connect_work"]])
        
        ## Get names back
        rownames(iol[["connect_total"]]) <- rows
        colnames(iol[["connect_total"]]) <- cols
      
        rownames(iol[["connect_work"]]) <- rows
        colnames(iol[["connect_work"]]) <- cols
        
    }

    return(iol)

}

################################################################################
#' Initialization of reference data.
#'
#' The function adds a component to the optimization parameter list passed in
#' as parameter \code{op}. The added component \code{opt.target} contains
#' observed data depending which data are provided on input to the function
#' [load.input()]. The function additionally checks whether execution of the
#' optimization procedure is possible based on the selected optimization
#' targets and the provided data.
#'
#' @param iol Input data list. Use [load.input()] to load needed fies and
#'            [init.connectivity()] in order to create a valid date layout.
#' @param op List with steering parameters for the optimization process. \cr
#'           Use [set.optimization.params()] in order to create a valid layout.
#' @param sp List with static model parameters. Use [set.static.params()] to 
#'           create a valid layout.
#' @param sim.struc List with population data. Use [init.spatial.population()]
#'                  in order to create a valid layout.
#'
#' @return The list with steering parameters for the optimization process passed
#'         in as parameter \code{op} with an additional component
#'         \code{opt.target} carying observed data, prepared to be used as
#'         target data in the optimization procedure of the [CoSMic()]
#'         function.
#'
#' @section ToDo:
#' Implement ot\[\[dea.nuts2\]\]
#'
#' @export
init.reference.data  <- function (iol, op, sp, sim.struc) {

    ## List of optimization targets ----------------------------------
    ot <- list()

    ## Reported deaths #########################################################
    
    ot[["dea"]]     <- iol$dead.cases.by.country
    
    ot[["dea.bs"]]  <- iol$dead.cases.by.state
   
    ## If country wide deaths are not given but opt target is country-deaths ---
    if (op$opt.target.deaths &
        (op$opt.target.region == "country") &
        is.null(ot$dea)) {
        stop(paste("Deaths are selected as global optimization target but",
                   "no observerd data are given by dead.cases.by.country"))
    }
    
    if ( ! is.null(ot$dea) ) {
        
        ## Subset to relevant timeframe from seed_date to ----------------------
        ## end of list - 10 days. to account for reporting delay ---------------
        ot$dea <- ot$dea[ot$dea$date %in%
                         seq(from = as.Date(sp$seed_date)+1,
                             to   = max(ot$dea$date) - 10,
                             by   = "days"),
                         ]
    }
    
    ## If deaths by state are not given but opt target is state-deaths ---------
    if (op$opt.target.deaths &
        (op$opt.target.region == "state") &
        is.null(ot$dea.bs)) {
        stop(paste("Deaths are selected as local optimization target but",
                   "no observerd data are given by dead.cases.by.state"))
    }
    
    if ( ! is.null(ot$dea.bs) ) {
        
        ## Subset to relevant timeframe from seed_date to ----------------------
        ## end of list - 10 days. to account for reporting delay ---------------
        ot$dea.bs <- ot$dea.bs[ot$dea.bs$date %in%
                         seq(from = as.Date(sp$seed_date)+1,
                             to   = max(ot$dea.bs$date) - 10,
                             by   = "days"),
                         ]
    }

    ## Reported ICU cases ######################################################
    
    ot[["icu"]]        <- iol$icu.cases.by.country
    
    ot[["icu.bs"]]     <- iol$icu.cases.by.state

    ot[["icu.nuts2"]]  <- iol$icu.cases.by.county
    
    ## If country wide ICU cases are not given but opt target is country-ICU ---
    if (op$opt.target.icu &
        (op$opt.target.region == "country") &
        is.null(ot$icu)) {
        stop(paste("ICU cases are selected as global optimization target but",
                   "no observerd data are given by icu.cases.by.country"))
    }
    
    if ( ! is.null(ot$icu) ) {
        
        ## Subset to relevant timeframe from seed_date to end of list ----------
        ot$icu <- ot$icu[ot$icu$date %in%
                         seq(from = as.Date(sp$seed_date)+1,
                             to   = max(ot$icu$date),
                             by   = "days"),
                         ]
    }
    
    ## If icu by state is not given but opt target is state icu ----------------
    if (op$opt.target.icu &
        (op$opt.target.region == "state") &
        is.null(ot$icu.bs)) {
        stop(paste("ICU cases are selected as local optimization target but",
                   "no observerd data are given by icu.cases.by.state"))
    }
    
    if ( ! is.null(ot$icu.bs) ) {
        
        ## Subset to relevant timeframe from seed_date to end of list ----------
        ot$icu.bs <- ot$icu.bs[ot$icu.bs$date %in%
                         seq(from = as.Date(sp$seed_date)+1,
                             to   = max(ot$icu.bs$date),
                             by   = "days"),
                         ]
    }

    ## If icu by county is not given but opt target is nuts2 icu ----------------
    if (op$opt.target.icu &
        (op$opt.target.region == "nuts2") &
        is.null(ot$icu.nuts2)) {
        stop(paste("ICU cases are selected as local optimization target but",
                   "no observerd data are given by icu.cases.by.county"))
    }

    if ( ! is.null(ot$icu.nuts2) ) {
        
        ## Subset to relevant timeframe from seed_date to end of list ----------
        ot$icu.nuts2 <- ot$icu.nuts2[ot$icu.nuts2$date %in%
                                     seq(from = as.Date(sp$seed_date)+1,
                                         to   = max(ot$icu.nuts2$date),
                                         by   = "days"),
                                     ]

        ## Merge in Nuts2 column from iol$counties -----------------------------
        ot$icu.nuts2 <- merge(iol$counties[,c("dist_id","Nuts2")],ot$icu.nuts2,
                              by="dist_id")

        ot$icu.nuts2$Nuts2 <- trimws(ot$icu.nuts2$Nuts2,which="both")
        ## Group first by NUTS2 regions and then by date and sum up ------------
        ot$icu.nuts2 <- ot$icu.nuts2 %>% group_by(Nuts2,date) %>% summarize(cases=sum(cases))

    }
    
    op$opt.target <- ot
    return(op)

}
