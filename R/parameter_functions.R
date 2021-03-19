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
#' Setup of execution parameters
#'
#' @param exec.procedure Set the execution procedure. Valid values are
#'                       "Basic-Param" or "Optimization"\cr
#'                       *Defaults to:* \code{"Basic-Param"}.
#' @param  parallel.method Set the parallelization method. Valid values are
#'                         "OMP", "MPI" or "PSOCK"\cr
#'                         *Defaults to:* \code{"OMP"}.
#' @param max.cores Set the maximum number of cores used in case
#'                  parallel.method = "OMP".\cr
#'                  *Defaults to:* \code{4}
#' @param omp.cluster.dbg Whether std.out from workers should be captured to
#'                        a file called cl.out.\cr
#'                        *Defaults to:* \code{FALSE}
#' @param data.dir Path to the directory from which input files are read.\cr
#'                 *Defaults to:* \code{"data"}
#' @param model.version The model version string.\cr
#'                      *Currently defaults to:* \code{12.0}
#' @param export_name File name addition for output files.\cr
#'                    *Defaults to:* \code{<model.version>-<YYYY-MM-DD_hh:mm:ss>}
#' 
#' @return A list with parameters needed to set up the execution of the CoSMic
#'         function. The default structure is:
#' ```
#' $exec.procedure
#' [1] "Basic-Param"
#' $parallel.method
#' [1] "OMP"
#' $max.cores
#' [1] 4
#' $omp.cluster.dbg
#' [1] FALSE
#' $data.dir
#' [1] "Data"
#' $model.version
#' [1] "12.0"
#' $export_name
#' [1] "v12.0-2020-11-07_21:53:00"
#' ```
#' @import rlist
#' @export
set.exec.params <- function(exec.procedure  = "Basic-Param",
                            parallel.method = "OMP",
                            max.cores       = 4,
                            omp.cluster.dbg = FALSE,
                            data.dir        = "data",
                            output.dir      = NULL,
                            model.version   = "12.0",
                            export_name     = NULL,
                            cp.write        = FALSE,
                            cp.time         = 0,
                            cp.reload       = FALSE,
                            cp.dir          = NULL ) {
    ep <- list()
    
    ## E1. Execution procedure ---------------------------------------
    ##  Valid values : "Basic-Param" or "Optimization"
    if (exec.procedure == "Basic-Param" |
        exec.procedure == "Optimization" ) {
        
        ep <- list.append(ep,
                          exec.procedure  = exec.procedure)
    } else {
        stop(paste("Unsupported value for exec.procedure specified:",
                   exec.procedure,"\n ",
                   "Please specify either \"Basic-Param\" or \"Optimzation\""))
    }

    ## E2. Parallel Execution Method ---------------------------------
    ##  Valid values : "OMP", "MPI" or "PSOCK"
    if (parallel.method %in% c("OMP","MPI","PSOCK")) {

        ## Optimization can not be used with Rmpi --------------------
        if (exec.procedure == "Optimization" &
            parallel.method == "MPI" ) {
            stop(paste("parallel.method=\"MPI\" can not be used wih",
                       "exec.procedure=\"Optimzation\"","\n ",
                       "Please specify one of \"OMP\" or \"PSOCK\""))
        }
        
        ep <- list.append(ep,
                          parallel.method=parallel.method)
    } else {
        stop(paste("Unsupported value for parallel.method specified:",
                   parallel.method,"\n ",
                   "Please specify one of \"OMP\", \"MPI\" or \"PSOCK\""))
    }

    ## Treat max.cores -----------------------------------------------
    ep <- list.append(ep,
                      max.cores=max.cores)

    ## Treat omp.cluster.dbg -----------------------------------------
    if (typeof(omp.cluster.dbg) != "logical") {
        stop("omp.cluster.dbg can only be TRUE or FALSE")
    }
    ep <- list.append(ep,
                      omp.cluster.dbg =omp.cluster.dbg )

    ## Data directory ------------------------------------------------
    if ( !dir.exists(data.dir) ) {
        warning(paste("Data dir",data.dir,"does not exist!"))
    }
    ep <- list.append(ep,
                      data.dir=data.dir)

    ## Treat model.version -------------------------------------------
    ep <- list.append(ep,
                      model.version=model.version)

    ## Treat export_name ----------------------------------------------
    if (is.null(export_name)) {
        export_name <- paste0("v",model.version,"-",
                               format.POSIXct(Sys.time(), format = "%Y-%m-%d_%H:%M:%S"))
    }
    
    ep <- list.append(ep, export_name=export_name)

    ## Output directory ----------------------------------------------
    if (is.null(output.dir)) {
        output.dir <- paste0("Results-",ep$export_name,"/")
    } else {
        if (length(grep("/$",output.dir)) == 0) {
            output.dir <- paste0(output.dir,"/")
        }
    }
        
    if ( !dir.exists(output.dir) ) {
        dir.create(output.dir)
    }
    ep <- list.append(ep,
                      output.dir=output.dir)

    ## Read checkpoint -----------------------------------------------
    if (typeof(cp.write) != "logical") {
        stop("cp.write can only be TRUE or FALSE")
    }
    ep <- list.append(ep,
                      cp.write = cp.write )

    ## Read checkpoint -----------------------------------------------
    if (class(cp.time) != "numeric") {
        stop("cp.time can only be a numeric value")
    }
    ep <- list.append(ep,
                      cp.time = cp.time )
    
    ## Reload checkpoint -----------------------------------------------
    if (typeof(cp.reload) != "logical") {
        stop("cp.reload can only be TRUE or FALSE")
    }
    ep <- list.append(ep,
                      cp.reload = cp.reload )
    ## Checkpoint directory ------------------------------------------
    if (is.null(cp.dir) & cp.reload) {
        stop(paste("cp.reload is TRUE but no cp.dir was given.",
                   "This is not possible."))
    }
    ep <- list.append(ep,
                      cp.dir = cp.dir )
    
    return(ep)
}

################################################################################
#' Function to save the current list of execution parameters.
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' 
#' @export
save.exec.params <- function(ep) {
    saveRDS(file=paste(ep$output.dir,
                    "exec.params-",ep$export_name,".RDS",sep=""),
            object=ep)
}

################################################################################
#' Setup of parameters in parameter space
#'
#' The function adds an element to the parameter space list pspace
#' 
#' @param param The name of the parameter to be set.
#' @param values The values of the parameter to be set.
#' @param type The parameter type. Allowed values are \code{direct} or
#'             \code{dist}.\cr
#'             *Defaults to:* \code{direct}
#' @param s.dev Deviations of the values in case of parameter type
#'              \code{dist}.\cr
#'              *Defaults to:* \code{NULL}
#' 
#' @return The function operates on the global scope and modifies the 
#'         parameter list pspace.\
#'         
#' @export
set.pspace <- function(param, values, type="direct", s.dev=NULL) {

    ## If pspace does not exist set it in global scope ---------------
    if ( !exists("pspace") ) {
        pspace <<- list()
    }

    ## For param only type character is supported --------------------
    if (typeof(param) != "character") {
        stop(paste("Given value for parameter param is not of type \"character\""))
    }
    ## For param only a single name is supported ---------------------
    if (length(param) != 1) {
        stop(paste("Given value for parameter param is not of length 1\n ",
                   "Please specify only a single parameter name."))
    }

    ## Check whether a supported parameter type is given -------------
    if ( ! (type %in% c("direct", "dist") ) ) {
        stop(paste("Given value for parameter type is not supported","\n ",
                   "Please specify one of \"direct\" or \"dist\""))
    }

    ## In case s.dev is present set type to "dist" -----------------------------
    if ( !is.null(s.dev) ) type  <- "dist"
    
    ## In case value is a data frame and sd is not given -----------------------
    ## reset type to directv 
    if ( (class(values) == "data.frame") & (is.null(s.dev)) ) {
        print("values was recognised as data.frame with sd not given. type is reset to 'directv'")
        type = "directv"
    }

    ## In case value is a data frame and sd is given ---------------------------
    ## reset type to distv 
    if ( (class(values) == "data.frame") & (!is.null(s.dev)) ) {
        print("values was recognised as data.frame with sd given. type is reset to 'distv'")
        type = "distv"
    }
    
    ## In case value is a list reset type to directl -----------------
    if (class(values) == "list") {
        print("values was recognised as list. type is reset to 'directl'")
        type = "directl"
    }
    
    ## If param is already in pspace isse a warning and replace ------
    if (param %in% names(pspace)) {
        warning(paste("Parameter",param,"is already set in pspace.\n ",
                      "Its values will be replaced."))

        pspace[[param]] <<- list(param = values,
                                 type  = type   )

    } else {
    
        pspace <<- list.append(pspace,
                               list(param = values,
                                    type  = type   )
                               )
        le <- length(pspace)
        
        names(pspace)[le] <<- param
    }

    ## If parameter type is "dist" set standard deviation ----------------------
    if (type %in% c("dist","distv")) {
        
        ## if s.dev is not given, stop -------------------------------
        if (is.null(s.dev)){
            stop(paste("type=\"dist\" was specified without s.dev\n ",
                       "This case is not supported."))
        }

        if (length(unlist(s.dev)) != length(unlist(values))) {
            stop(paste("length(unlist(s.dev)) != length(unlist(values) for param",param))
        }

        pspace[[param]]$sd    <<- s.dev
    }
  
    ## In case a file list is given as values load the data -----------
    if ( (type == "directl") &
         ( all(lapply(values,class) == "character") ) ) {
        print("List of character elements detected. Data will be loaded from file")
        jj  <- 1
        for ( ii in values ) {
            pspace[[param]]$param[[jj]] <<- read.table(file=ii)
            jj <- jj + 1
        }
    }
}
##
## The pspace list
## ================
##
## Inputs for parameter studies are stored in the list pspace. Prototypes of
## supported parameter notations in the list pspace are explained below.
##
## -----------------------------------------------------------------------------
## Type: "dist"
## ------------
##       This type is currently only working for floating point parameters,
##       since randomLHS returns floating point values and scaling to integers
##       is not done at the moment. It has to be speciifed as follows:
##
##       pspace <- list.append(pspace,
##                              <parameter_name>=list(
##                              param    = <value>,
##                              type     = "dist",
##                              sd       = <varition_in_percent>,
##                              num.type = <numerical_type> ))
##
##       The first line is equal for all parameters.
##       <parameter_name>       : Name of the parameter.
##       <value>                : Mean value or vector of mean values for vector
##                                parameters, around which to vary.
##       <variation_in_percent> : Max. variation or vector of max. variations of
##                                the parameter in percent of <value> i.e. the
##                                corresponding vector element of value. Set
##                                <variation_in_percent>=0 for each element that
##                                should not be varied.
##       <numerical_type>       : Currently only "num" is supported
##
## Examples:
## ---------
## 1.) Scalar Parameter R0: The mean Value 2.5 is varied max. 5% in positive
##                          and negative direction.
##
##     pspace <- list.append(pspace,
##                           R0=list(
##                           param    = 2.5,
##                           type     = "dist",
##                           sd       = 5,
##                           num.type = "num"))
##
## 2.) Vector Parameter R0effect: The first three elements of the param vector
##                                are varied by max .5% each in positive and 
##                                negative direction. The next two are varied
##                                by 2.5% and the sixth parameter is kept
##                                constant since the coresponding element in
##                                sd is equal 0.
##
##     pspace <- list.append(pspace,
##                           R0effect = list(
##                           param    = c(0.65,0.65,0.501,0.545,0.567,0.556),
##                           type     = "dist",
##                           sd       = c(   5,   5,   5 , 2.5 , 2.5 ,   0 ),
##                           num.type = "num"))
##
## -----------------------------------------------------------------------------
## Type: "direct"
## --------------
##       This type has to be speciifed as follows:
##
##       pspace <- list.append(pspace,
##                              <parameter_name>=list(
##                              param    = <value>,
##                              type     = "direct"))
##
##       The first line is equal for all parameters.
##       <parameter_name>       : Name of the parameter.
##       <value>                : Mean value or vector of values that
##                                should be evaluated.
## Examples:
## ---------
## 1.) Scalar Parameter icu_dur: The Prameter is evaluated with the values
##                               12 and a second time 14.
##     pspace <- list.append(pspace,
##                           icu_dur = list(
##                           param   =  c(12,14),
##                           type    = "direct"))
##
## 2.) Scalar Parameter icu_dur: The Prameter is evaluated with value 14
##
##     pspace <- list.append(pspace,
##                           icu_dur = list(
##                           param   =  14,
##                           type    = "direct"))
##
## -----------------------------------------------------------------------------
## Type: "directv"
## ---------------
##       This type has to be speciifed as follows:
##
##       pspace <- list.append(pspace,
##                             <parameter_name> = list(
##                             param    = data.frame(<value_columns>),
##                             type     = "directv"))
##
##       The first line is equal for all parameters.
##       <parameter_name>       : Name of the parameter.
##       <value_columns>        : data.frame columns. One for each vector
##                                element and one row for each value that 
##                                should be evaluated.
## Examples:
## ---------
## 1.) Vector parameter R0effect: 
##
##     pspace <- list.append(pspace,
##                      R0effect = list(
##                      param    = data.frame(R0effect1=0.65,
##                                            R0effect2=0.65,
##                                            R0effect3=0.501,
##                                            R0effect4=0.545,
##                                            R0effect5=0.567,
##                                            R0effect6=seq(0.4,1.,0.3 )),
##                      type     = "directv"))
##
##     R0effect1-5 are kept constant while R0effect takes the values
##     0.4, 0.7 and 1.0.
##
## -----------------------------------------------------------------------------
##
## The lhc variable
## ----------------
##
## In the section ### Prepare parameter space ### directly before the parallel
## foreach loop, the parameters stored in the pspace list are premutated against
## each other and the iter parameter and transferred to a data.frame whose rows
## each caries a set of parameters with which the model will be evaluated.
##
## lhc is printed right befor the execution of the foreach loop so that the user
## has an overview about the parameterspace to be evaluated.
##
## -----------------------------------------------------------------------------
##
## Implementation strategie for new scalar parameters:
## ---------------------------------------------------
##
## To replace a scalar parameter within the model code the original parameter
## has to be replaced by its named column in the lhc variable carying all
## parameters in pspace.
##
## <param_name> :  lhc[it.ss,"<param_name>"]
##
## E.g: sam_size has to be replaced with lhc[it.ss,"sam_size"]
##
## Implementation strategie for new vector parameters:
## ---------------------------------------------------
## To replace a vector parameter within the model code the original parameter
## has to be replaced by its named columns in the lhc variable carying all
## parameters in pspace.
##
## <param_name> :  lhc[it.ss,paste0("<_param_name>", <param_no>)]
##
## E.g: R0effect has to be replaced with lhc[1,paste0("R0effect",change)]
##

################################################################################
#' Function to save the current psapce list. 
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' @param pspace The parameter list pspace.
#' 
#' @export
save.pspace <- function(ep, pspace) {
    saveRDS(file=paste(ep$output.dir,
                    "pspace-",ep$export_name,".RDS",sep=""),
            object=pspace)
}

################################################################################
#' Setup of static parameters
#'
#' @export
set.static.params <- function(pspace,
                              seed.in.inner.loop  = FALSE,
                              seed.base           = NULL,             
                              country             = "Germany",
                              restrict            = TRUE,
                              sim.regions         = c("Schleswig-Holstein","Hamburg","Niedersachsen","Bremen"),
                              sam_prop.ps         = c(1.0,1.0,1.0,1.0),
                              sim_pop             = "proportional",
                              ini_infected        = 10,
                              seed_infections     = "data",
                              seed_date           = "2020-03-09",
                              seed_before         = 7,
                              time_n              = NULL,
                              inf_dur             = 3,
                              cont_dur            = 2,
                              ill_dur             = 8,
                              icu_per_day         = c(0,0,0,0,0,0,0,8),
                              less_contagious     = 0.7,
                              R0_force            = 0,
                              immune_stop         = TRUE,
                              import_R0_matrix    = FALSE,
                              R0change            = lapply(seq(1,by=7,length.out=20),
                                                           function(x){c(x,x+6)}),
                              R0county            =  as.list(rep("ALL",20)),
                              R0delay             = TRUE,
                              R0delay_days        = 5,
                              R0delay_type        = "linear",
                              endogenous_lockdown = FALSE,
                              lockdown_effect     = 0.39,
                              lockdown_connect    = 0.5,
                              lockdown_threshold  = 100,
                              lockdown_days       = 10,
                              control_age_sex     = "age",
                              iter                = 4,
                              lhc.samples         = NULL,
                              lhc.reload          = FALSE,
                              gplots              = FALSE,
                              cplots              = FALSE,
                              cplots.states       = FALSE,
                              cplots.nuts2        = FALSE,
                              results             = "Reduced",
                              sp.states           = NULL) {
    
    sp <- list()

    ## Treat seed.in.inner.loop --------------------------------------
    if (typeof(seed.in.inner.loop) != "logical") {
        stop("seed.in.inner.loop can only be TRUE or FALSE")
    }
    sp <- list.append(sp,
                      seed.in.inner.loop=seed.in.inner.loop)

    ## Treat seed.base -----------------------------------------------
    if ( seed.in.inner.loop & is.null(seed.base) ) {
        print(paste("seed.in.inner.loop is TRUE with seed.base is NULL.",
                    "This is not possible. seed.base is set to 42."))
        seed.base <- 42
    }
    sp <- list.append(sp,
                      seed.base=seed.base)
    
    ## Treat country -------------------------------------------------
    if ( ! (country %in% "Germany") ) {
        stop("Currently only Germany is supported")
    }
    sp <- list.append(sp,
                      country=country)

    ## Treat restrict ------------------------------------------------
    if (typeof(restrict) != "logical") {
        stop("restrict can only be TRUE or FALSE")
    }
    sp <- list.append(sp,
                      restrict=restrict)
    
    ## Treat sim.regions ---------------------------------------------
    if ( ! restrict & ! is.null( sim.regions )) {
        warning(paste("sim.regions = ",paste(sim.regions,collapse=","),"\n",
                      "this will have no effect since restrict = FALSE.","\n",
                      "This warning will disappear if sim.regions",
                      "is set to NULL"))
    }
    
    sp <- list.append(sp,
                 sim.regions=sim.regions)

    ## Treat sam_prop.ps ---------------------------------------------
    if ( any( sam_prop.ps <= 0.0 ) ) {
        stop(paste("sam_prop.ps contains elements which are <= 0",
                   "This is not possible!","\n",
                   "sam_prop.ps = [",paste(sam_prop.ps,collapse=","),"]"))
    }
    sp <- list.append(sp,
                      sam_prop.ps=sam_prop.ps)

    ## Treat sim_pop -------------------------------------------------
    if ( ! ( sim_pop %in% c("proportional", "random") ) ) {
        stop(paste("sim_pop = ",sim_pop,". This is not supported!",
                   "please select 'proportional' or 'random'"))
    }
    sp <- list.append(sp,
                      sim_pop = sim_pop)

    ## Treat ini_infected ---------------------------------------------
    if ( ini_infected < 0 ) {
        stop("ini_infected < 0 . This is not possible!")
    }
    sp <- list.append(sp,
                      ini_infected=ini_infected)

    ## Treat seed_infections -----------------------------------------
    if ( ! ( seed_infections %in% c("random", "county", "data") ) ) {
        stop(paste("seed_infections = ",seed_infections,". This is not supported!",
                   "please select 'random', 'county' or 'data'"))
    }
    sp <- list.append(sp,
                      seed_infections=seed_infections)

    ## Treat seed_date -----------------------------------------------
    sp <- list.append(sp,
                      seed_date=as.Date(seed_date))

    ## Treat seed_before -----------------------------------------------
    if ( seed_before < 0 ) {
        stop("seed_before < 0 . This is not possible!")
    }
    sp <- list.append(sp,
                      seed_before=seed_before)
    
    ## Treat inf_dur -------------------------------------------------
    if ( inf_dur < 0 ) {
        stop("inf_dur < 0 . This is not possible!")
    }
    sp <- list.append(sp,
                      inf_dur=inf_dur)

    ## Treat cont_dur ------------------------------------------------
    if ( cont_dur < 0 ) {
        stop("cont_dur < 0 . This is not possible!")
    }
    sp <- list.append(sp,
                      cont_dur=cont_dur)

    ## Treat ill_dur -------------------------------------------------
    if ( ill_dur < 0 ) {
        stop("ill_dur < 0 . This is not possible!")
    }
    sp <- list.append(sp,
                      ill_dur=ill_dur)
    
    ## Treat icu_per_day ---------------------------------------------
    if (length(icu_per_day)!=ill_dur)  {
        stop(paste("Length of icu_per_day = ",length(icu_per_day)," while",
                   "ill_dur =",ill_dur,". \n",
                   "length(icu_per_day)==ill_dur has to apply"))
    }
    if (mean(icu_per_day)!=1) {
        stop(paste("mean(icu_per_day) =",mean(icu_per_day),
                   "but it has to be equal to 1"))
    }
    sp <- list.append(sp,
                      icu_per_day=icu_per_day)

    ## Treat less_contagious -----------------------------------------
    if ( less_contagious < 0 ) {
        stop("less_contagious < 0 . This is not possible!")
    }
    if ( less_contagious > 1.0 ) {
        stop("less_contagious > 1.0 . This is not possible!")
    }
    sp <- list.append(sp,
                      less_contagious=less_contagious)

    ## Treat R0_force ------------------------------------------------
    if ( R0_force < 0 ) {
        stop("R0_force < 0 . This is not possible!")
    }
    if ( R0_force > 1.0 ) {
        stop("R0_force > 1.0 . This is not possible!")
    }
    sp <- list.append(sp,
                      R0_force=R0_force)

    ## Treat immune_stop ----------------------------------------------
    if (typeof(immune_stop) != "logical") {
        stop("immune_stop can only be TRUE or FALSE")
    }
    sp <- list.append(sp,
                      immune_stop=immune_stop)

    ## Treat import_R0_matrix ----------------------------------------
    sp <- list.append(sp,
                      import_R0_matrix=import_R0_matrix)

    ## Treat R0change ------------------------------------------------
    ## If R0_effect is already in pspace check again its length
    if ( "R0_effect" %in% names(pspace) ) {
        if ( class( pspace[["R0_effect"]]$param ) == "numeric" ) {

            len.R0_effect <- length(pspace[["R0_effect"]]$param)
            
        } else if ( class( pspace[["R0_effect"]]$param ) == "data.frame" ) {

            len.R0_effect <- dim(pspace[["R0_effect"]]$param)[2]
            
        } else if ( class( pspace[["R0_effect"]]$param ) == "list" ) {

            len.R0_effect <- dim(pspace[["R0_effect"]]$param[[1]])[1]

        }
        if (length( R0change ) > len.R0_effect) {
            stop("R0change is of different length than R0_effect.",
                 "This is not possible!")
        }
    }
    sp <- list.append(sp,
                      R0change=R0change)

    ## Treat time_n ---------------------------------------------------
    sp <- list.append(sp,
                      time_n=time_n)
    
    if (is.null(sp$time_n)) {
        sp$time_n <- max(unlist(sp$R0change))+1
    } else if (max(unlist(sp$R0change)) > sp$time_n) {
        sp$time_n <- max(unlist(sp$R0change))+1
    }
    
    ## Treat R0county ------------------------------------------------
    if (length(R0county) != length(R0change)) {
        stop("R0change is of different length than R0county.",
             "This is not possible!")
    }
    sp <- list.append(sp,
                      R0county=R0county)

    
    ## Treat R0delay -------------------------------------------------
    if (typeof(R0delay) != "logical") {
        stop("R0delay can only be TRUE or FALSE")
    }
    sp <- list.append(sp,
                      R0delay=R0delay)

    ## Treat R0delay_days --------------------------------------------
    if (typeof(R0delay_days) < 0) {
        stop("R0delay_days is < 0. This is not possible!")
    }
    sp <- list.append(sp,
                      R0delay_days=R0delay_days)

    ## Treat R0delay_type --------------------------------------------
    if ( ! (R0delay_type %in% c("logistic","linear") ) ) {
        stop("R0delay_type is neither 'logistic' nor 'linear'. This is not supported!")
    }
    sp <- list.append(sp,
                      R0delay_type=R0delay_type)  

    ## Treat endogenous_lockdown -------------------------------------
    if (typeof(endogenous_lockdown) != "logical") {
        stop("endogenous_lockdown can only be TRUE or FALSE")
    }
    sp <- list.append(sp,
                      endogenous_lockdown=endogenous_lockdown)

    ## Treat lockdown_effect -----------------------------------------
    if ( lockdown_effect < 0 ) {
        stop("lockdown_effect < 0 . This is not possible!")
    }
    if ( lockdown_effect > 1.0 ) {
        stop("lockdown_effect > 1.0 . This is not possible!")
    }
    sp <- list.append(sp,
                      lockdown_effect=lockdown_effect)
    
    ## Treat lockdown_connect -----------------------------------------
    if ( lockdown_connect < 0 ) {
        stop("lockdown_connect < 0 . This is not possible!")
    }
    if ( lockdown_connect > 1.0 ) {
        stop("lockdown_connect > 1.0 . This is not possible!")
    }
    sp <- list.append(sp,
                      lockdown_connect=lockdown_connect)
    
    ## Treat lockdown_threshold ---------------------------------------
    if ( lockdown_threshold < 0 ) {
        stop("lockdown_threshold < 0 . This is not possible!")
    }
    sp <- list.append(sp,
                      lockdown_threshold=lockdown_threshold)

    ## Treat lockdown_days -------------------------------------------
    if ( lockdown_days < 0 ) {
        stop("lockdown_days < 0 . This is not possible!")
    }
    sp <- list.append(sp,
                      lockdown_days=lockdown_days)

    ## Treat control_age_sex -----------------------------------------
    if ( !( control_age_sex %in% c("NONE", "age", "age_sex") ) ) {
        stop(paste("control_age_sex = ", contorl_age_sex,". This is not supported!\n",
                   "Please select 'NONE', 'age' or 'age_sex'"))
    }
    sp <- list.append(sp, control_age_sex=control_age_sex)
    
    ## Treat iter ----------------------------------------------------
    if (is.null(iter)) iter  <- 10
    sp <- list.append(sp, iter=iter)
    
    ## Treat lhc.samples ---------------------------------------------
    if (is.null(lhc.samples)) lhc.samples <- 1

    ## Extract parameter types from paramter list ---------------------
    pspace.types <- lapply(pspace,function(x){x$type})
    
    ## Extract standard deviations of type dist -----------------------
    pspace.sd <- lapply(pspace[which(pspace.types %in% c("dist","distv"))],
                        function(x){x$sd})

    if ((lhc.samples <= 1) & length(pspace.sd) > 0) {         
        stop(paste("Found parameters of dist type in pspace but lhc.samples <= 1",
                   "This is not meaningfull."))
    }
    
    sp <- list.append(sp, lhc.samples=lhc.samples)

    ## Treat lhc.reload ----------------------------------------------
    if (typeof(lhc.reload) != "logical") {
        stop("lhc.reload can only be TRUE or FALSE")
    }
    sp <- list.append(sp, lhc.reload=lhc.reload)
 
    ## Treat gplots --------------------------------------------------
    if (typeof(gplots) != "logical") {
        stop("gplots can only be TRUE or FALSE")
    }
    sp <- list.append(sp, gplots=gplots)

    ## Treat cplots --------------------------------------------------
    if (typeof(cplots) != "logical") {
        stop("cplots can only be TRUE or FALSE")
    }
    sp <- list.append(sp, cplots=cplots)

    ## Treat cplots.states -------------------------------------------
    if (typeof(cplots.states) != "logical") {
        stop("cplots.states can only be TRUE or FALSE")
    }
    sp <- list.append(sp, cplots.states=cplots.states)

    ## Treat cplots.nuts2 --------------------------------------------
    if (typeof(cplots.nuts2) != "logical") {
        stop("cplots.nuts2 can only be TRUE or FALSE")
    }
    sp <- list.append(sp, cplots.nuts2=cplots.nuts2)

    ## Treat results -------------------------------------------------
    sp <- list.append(sp, results=results)
    if ((sp$results != "ALL") &
        !all(!sp$gplots,!sp$cplots,!sp$cplots.states,sp$cplots.nuts2) ) {

        ## Treat gplots --------------------------------------------------
        sp$gplots <- FALSE
        
        ## Treat cplots --------------------------------------------------
        sp$cplots <- FALSE
        
        warning(paste("sp$results = ",sp$results," plotting on country level was disabled."))
        
    }

    ## Treat speaking names of health states -----------------------------------
    if (is.null(sp.states) | (length(sp.states) == 7)) {

        if ((length(sp.states) != 7) & (!is.null(sp.states))) {
            warning("The CoSMic SIER model operates with 7 health conditions.",
                    "Please provide as many speaking condition names")
        }
        
        sp.states  <- c("Healthy / Never infected",
                        "Infected, not contagious",
                        "Infected, contagious",
                        "Ill, contagious",
                        "Ill, ICU",
                        "Immune",
                        "Dead")
        
        sp <- list.append(sp, sp.states=sp.states)
        
    } else {

        sp <- list.append(sp, sp.states=sp.states)
        
    } 

    return(sp)
}

################################################################################
#' Function to save the current list of satic model parameters. 
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' @param sp A list with static model parameters as described in [set.static.params()].
#' 
#' @export
save.static.params <- function(ep, sp) {
    saveRDS(file=paste(ep$output.dir,
                    "static.params-",ep$export_name,".RDS",sep=""),
            object=sp)
}

################################################################################
#'  Loading input data
#'
#' @export
load.input <- function(data.dir              = "./",
                       trans.pr              = NULL,
                       pop.data              = NULL,
                       inf.cases             = NULL,
                       dead.cases            = NULL,
                       connect.total         = NULL,
                       connect.work          = NULL,
                       sts                   = NULL,
                       cnts                  = NULL,
                       R0.matrix.inp         = NULL,
                       dead.cases.by.state   = NULL,
                       dead.cases.by.country = NULL,
                       icu.cases.by.county   = NULL,
                       icu.cases.by.state    = NULL,
                       icu.cases.by.country  = NULL,
                       lhc.data              = NULL	) {
    
    iol <- list()

    ## Transition probabilities ------------------------------------------------
    if (is.null(trans.pr)) {
        data(trans_pr)
        iol[["trans_pr"]]  <- trans_pr
    } else {
        iol[["trans_pr"]] <- read.csv(paste(data.dir,trans.pr,sep="/"))
    }
    
    ## Population structure ----------------------------------------------------
    if (is.null(pop.data)) {
        data("181231_pop_age_sex_distr_ger")
        iol[["pop"]]  <- pop
    } else {
        iol[["pop"]]  <- read.csv(paste(data.dir,pop.data,sep="/"))
    }
    
    ## Seed data for infection seeding during model startup --------------------
    if (is.null(inf.cases)) {
        data(infections)
        iol[["seed"]]  <- seed
    } else {
        iol[["seed"]] <- read.csv(paste(data.dir,inf.cases,sep="/"))
        }

    ## Seed data for seedinf of dead cases during model startup ----------------
    if (is.null(dead.cases)) {
        data(dead_cases)
        iol[["seed_dea"]]  <- seed_dea
    } else {
        iol[["seed_dea"]] <- read.csv(paste(data.dir,dead.cases,sep="/"))
    }
    
    ## Connectivity matrix: Total population -----------------------------------
    if (is.null(connect.total)) {
        data(connect_total)
        iol[["connect_total"]]  <- connect_total
    } else {
        iol[["connect_total"]] = read.csv(paste(data.dir,connect.total,sep="/"))
    }
    
    ## Connectivity matrix: Working population ---------------------------------
    if (is.null(connect.work)) {
        data(connect_work)
        iol[["connect_work"]]  <- connect_work
    } else {
        iol[["connect_work"]] = read.csv(paste(data.dir,connect.work,sep="/"))
    }

    ## State information -------------------------------------------------------
    if (is.null(sts)) {
        data(states)
        iol[["states"]]  <- states
    } else {
        iol[["states"]] <- read.csv(paste(data.dir,sts,sep="/"),
                           stringsAsFactor=FALSE)
    }
    
    ## County information ------------------------------------------------------
    if (is.null(cnts)) {
        data(counties)
        iol[["counties"]]  <- counties
    } else {
        iol[["counties"]] <- read.csv(paste(data.dir,cnts,sep="/"),
                             stringsAsFactor=FALSE)
    }

    ## Change class of date column of seed to class Date -----------------------
    if ( ! ("date" %in% names(iol$seed)) ) {
        stop(paste("Colname \"date\" is missing in file \n",
                   inf.cases))
    }
    iol$seed$date <- as.Date(iol$seed$date)

    ## Load R0 data as matrix  -----------------------------------------------------------
    if ( !is.null(R0.matrix.inp) ) {
        iol <- list.append(iol,
                           R0_raw = read.csv(paste(data.dir,
                                                   R0.matrix.inp, sep=""))
                           )
    }

    ## Load reference data for dead cases for the complete country -------------
    if ( ! is.null( dead.cases.by.country ) ) {
        
        iol <- list.append(
            iol, dead.cases.by.country =
                     read.csv(file = paste(data.dir,
                                           dead.cases.by.country,sep="/")))
        
        ## Check whether the relevant columns exist ----------------------------
        if ( ! all(c("date","cases") %in% names(iol$dead.cases.by.country) ) ) {
            stop("Please check the structure of the file given by dead.cases.by.country \n ",
                 "Column names found are :",
                 paste(names(iol$dead.cases.by.country),collapse=","),"\n ",
                 "But we need 'date' and 'cases'.")
        }
        
        iol$dead.cases.by.country$date <- as.Date(iol$dead.cases.by.country$date)
        
    } else {
        iol <- list.append(iol, dead.cases.by.country = NULL)
    }
    
    ## Load reference data for dead cases by state -----------------------------
    if ( ! is.null( dead.cases.by.state ) ) {
        
        iol <- list.append(
            iol, dead.cases.by.state =
                     read.csv(file = paste(data.dir,
                                           dead.cases.by.state,sep="/")))
        
        ## Check whether the relevant columns exist ----------------------------
        if ( ! all(c("date","SummeTodesfaelle","Bundesland") %in% names(iol$dead.cases.by.state) ) ) {
            stop("Please check the structure of the file given by dead.cases.by.state \n ",
                 "Column names found are :",
                 paste(names(iol$dead.cases.by.state),collapse=","),"\n ",
                 "But we need 'date','SummeTodesfaelle', 'Bundesland' .")
        }
        
        iol$dead.cases.by.state$date <- as.Date(iol$dead.cases.by.state$date)
        ## Rename Columns to English names -------------------------------------
        names(iol$dead.cases.by.state)[
            which(names(iol$dead.cases.by.state)=="SummeTodesfaelle")] <- "deaths"
        names(iol$dead.cases.by.state)[
            which(names(iol$dead.cases.by.state)=="Bundesland")] <- "state"
        
    } else {
        iol <- list.append(iol, dead.cases.by.state = NULL)
    }
    
    ## Load reference data for ICU cases for the complete country --------------
    if ( ! is.null( icu.cases.by.country ) ) {

        iol <- list.append(
            iol, icu.cases.by.country =
                     read.csv(file = paste(data.dir,
                                           icu.cases.by.country,sep="/")))

        ## Check whether the relevant columns exist ----------------------------
        if ( ! all(c("date","cases") %in% names(iol$icu.cases.by.country) ) ) {
            stop("Please check the structure of the file given by icu.cases.by.country \n ",
                 "Column names found are :",
                 paste(names(iol$icu.cases.by.country),collapse=","),"\n ",
                 "But we need 'date' and 'cases'.")
        }

        iol$icu.cases.by.country$date <- as.Date(iol$icu.cases.by.country$date)
        
    } else {
        iol <- list.append(iol, icu.cases.by.country = NULL)
    }
    
    ## Load reference data for icu cases by state ------------------------------
    if ( ! is.null( icu.cases.by.state ) ) {

        iol <- list.append(
            iol, icu.cases.by.state =
                     read.csv(file = paste(data.dir,
                                           icu.cases.by.state,sep="/")))

        ## Check whether the relevant columns exist ----------------------------
        if ( ! all(c("date","cases","state") %in% names(iol$icu.cases.by.state) ) ) {
            stop("Please check the structure of the file given by icu.cases.by.state \n ",
                 "Column names found are :",
                 paste(names(iol$icu.cases.by.state),collapse=","),"\n ",
                 "But we need 'date' and 'cases' and 'state'.")
        }

        iol$icu.cases.by.state$date <- as.Date(iol$icu.cases.by.state$date)

    } else {
        iol <- list.append(iol, icu.cases.by.state = NULL)
    }

    ## Load reference data for icu cases by county -----------------------------
    if ( ! is.null( icu.cases.by.county ) ) {

        iol <- list.append(
            iol, icu.cases.by.county =
                     read.csv(file = paste(data.dir,
                                           icu.cases.by.county,sep="/")))

        ## Check whether the relevant columns exist ----------------------------
        if ( ! all(c("date","cases","dist_id") %in% names(iol$icu.cases.by.county) ) ) {
            stop("Please check the structure of the file given by icu.cases.by.county \n ",
                 "Column names found are :",
                 paste(names(iol$icu.cases.by.county),collapse=","),"\n ",
                 "But we need at least 'date' and 'cases' and 'dist_id'.")
        }

        
        iol$icu.cases.by.county$POSIXct <- as.POSIXct(iol$icu.cases.by.county$date)
        iol$icu.cases.by.county$date    <- as.Date(iol$icu.cases.by.county$date)

    } else {
        iol <- list.append(iol, icu.cases.by.county = NULL)
    }

    ## Load prepared latin hypercube data --------------------------------------
    if ( ! is.null( lhc.data ) ) {
        iol <- list.append(
            iol, lhc.data = read.table(file = paste(data.dir,
                                                    lhc.data,sep="/")))
    } else {
        iol <- list.append(iol, lhc.data = NULL)
    }
    
    return(iol)
}

################################################################################
#' Function to save the current list of loaded input data.
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' @param iol A list with loaded input data as described in [load.input()].
#' 
#' @export
save.input <- function(ep, iol) {
    saveRDS(file=paste(ep$output.dir,
                    "input-",ep$export_name,".RDS",sep=""),
            object=iol)
}

################################################################################
#'  Setup of optimization parameters
#'
#' @export
set.optimization.params <- function(opt.target.icu    ,  
                                    opt.target.deaths ,   
                                    opt.target.region ,  
                                    opt.names         ,  
                                    opt.lb            ,  
                                    opt.ub            ,
                                    opt.pop.size      ,
                                    opt.max.iter      ,  
                                    use.sug.sol       ,
                                    ep, sp, pspc, opt.filter=NULL ) {

    ## In case optimization is requested at least one target has to be selected
    if ( ep$exec.procedure == "Optimization" ) {
        if ( ! ( opt.target.icu | opt.target.deaths ) ) {
            stop(paste("Neither opt.target.icu nor opt.target.deaths are ",
                       "TRUE with ep$exec.procedure == 'Optimization' \n",
                       "This will not work !"))
        }

        lhc <- init.lhc(pspc, sp)
        if ( dim(lhc)[1] > sp$iter ) {
            stop(paste0("It appears that some parameters in pspace are set to",
                        "be variable since dim(lhc)[1] > sp$iter. This is not",
                        "possible with ep$exec.procedure = 'Optimization' !"))
        }
        
    }

    ## Check opt.target.region for valid values --------------------------------
    if ( !(opt.target.region %in% c("country","state","nuts2") ) ) {
        stop(paste("opt.target.region is ",opt.target.region,
                   "this value is not supported. Please select one of ",
                   "'country','state' or,'nuts2'"))
    }

    ## Check whether opt.names, opt.lb and opt.ub are of eual length -----------
    if ( length(opt.lb) != length(opt.names)  ) {
        stop(paste("length(opt.lb) != length(opt.names) \n",
                   "This is not possible"))
    }
    if ( length(opt.ub) != length(opt.names)  ) {
        stop(paste("length(opt.ub) != length(opt.names) \n",
                   "This is not possible"))
    }
    if ( length(opt.ub) != length(opt.lb) ) {
        stop(paste("length(opt.ub) != length(opt.lb) \n",
                   "This is not possible"))
    }

    ## Treat opt.filter --------------------------------------------------------
    if (is.null(opt.filter)) {
        opt.filter <- TRUE
    }
    if (typeof(opt.filter) != "logical") {
        stop("opt.filter can only be TRUE or FALSE")
    }
    
    op <- list(
        opt.target.icu      = opt.target.icu,    
        opt.target.deaths   = opt.target.deaths ,
        opt.target.region   = opt.target.region ,
        opt.names           = opt.names         ,
        opt.lb              = opt.lb            ,
        opt.ub              = opt.ub            ,
        opt.pop.size        = opt.pop.size      ,
        opt.max.iter        = opt.max.iter      ,
        use.sug.sol         = use.sug.sol       ,
        opt.filter          = opt.filter)

    return(op)
}

################################################################################
#' Function to save the current list of loaded input data.
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' @param op A list with parameters steering the optimization procedure as
#'           described in [set.optimization.params()].
#' 
#' @export
save.optimization.params <- function(ep, op) {
    saveRDS(file=paste(ep$output.dir,
                    "opt.params-",ep$export_name,".RDS",sep=""),
            object=op)
}

################################################################################
#' Reload a checkpoint.
#'
#' The function loads data necessary to do a checkpoint restart and checks
#' them for usability and differences.
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' @param sp A list with static model parameters as described in [set.static.params()].
#' 
#' @export
checkpoint.check.reload <- function(ep, sp) {

    sp.file <- dir(ep$cp.dir)[which(grepl(pattern="static.params-(.*)RDS",dir(ep$cp.dir)))]
    if (length(sp.file) == 0) {
        stop(paste("File with static parameter list was not found during",
                    "Checkpoint reload check!","\n", "Please provide a file with ",
                    "naming pattern 'static.params-(.*)RDS' in ep$cp.dir\n",
                    ep$cp.dir))
    }
    if (length(sp.file) > 1) {
        stop(paste("Found more than one file with static parameter list ",
                   "during checkpoint reload check!\n", "Please provide a",
                   "directory for ep$cp.dir in which only one file with ",
                    "naming pattern 'static.params-(.*)RDS' exsists"))
    }
    cp.sp   <- readRDS(paste(ep$cp.dir,sp.file,sep="/"))

    if ( !all(names(sp[!(sp %in% cp.sp)]) %in%
              c("R0change","R0county","time_n","iter","lhc.samples","results")) ) {
        stop(paste("In checkpoint reload changes were detected in the following",
                    "elements of the static model parameter list.\n",
                    paste(names(sp[!(sp %in% cp.sp)]),collapse=","),"\n",
                    "This is not possible! Only changes in 'R0change','R0county',",
                    "'iter' ,'lhc.samples','results' and 'time_n' are allowed."))
    }
   
}

################################################################################
#' Map R0effects from NUTS-1 to NUTS-2
#'
#' The function maps R0effects on NUTS-1 i.e. German state level to R0effects
#' on NUTS-2 level
#'
#' @param ep An execution parameter list as decribed in [set.exec.params()].
#' @param sp A list with static model parameters as described in [set.static.params()].
#' 
#' @export
map.R0effects <- function(R0effect.nuts2,R0effect.states,rows=NULL) {

    if (is.null(rows)) {
        rows <- dim(R0effect.nuts2)[1]
    }

    if (dim(R0effect.nuts2)[1] < max(rows)) {
        stop(paste("dim(R0effect.nuts2)[1] < max(rows) this is not possible."))
    }
    if (dim(R0effect.states)[1] < max(rows)) {
        stop(paste("dim(R0effect.states)[1] < max(rows) this is not possible."))
    }
    
    map.names <- iol$states[unique(data.frame(
                         state_id=as.integer(iol$counties$dist_id/1000),
                         Nuts2=iol$counties$Nuts2))[,"state_id"],"Shortcut"]
    
    R0effect.nuts2[rows,] <- R0effect.states[rows,map.names]
    
}
