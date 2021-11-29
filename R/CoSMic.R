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
#' Function executing the simulation model.
#'
#' @param ep Execution parameter list. Use [set.exec.params()] in order to
#'           create a valid layout.
#' @param sp List with static model parameters. Use [set.static.params()] to 
#'           create a valid layout.
#' @param iol Input data list. Use [load.input()] to load needed fies and
#'            [init.connectivity()] in order to create a valid date layout.
#' @param pspace List holding the parameter space with potentially variable
#'               model parameters. Use the setter function [set.pspace()] to 
#'               add parameters.
#' @param sim.struc List with population data. Use [init.spatial.population()]
#'                  in order to create a valid layout.
#' @param op List with steering parameters for the optimization process. \cr
#'           Use [set.optimization.params()] in order to create a valid layout 
#'           and\cr [init.reference.data()] in order to init the optimization 
#'           targets based on observed data.
#' @param opt Numeric vector with model parameters subject to optimization.
#'
#' @return Depends upon the selected execution procedure given by
#'         ep$exec.procedure.
#'         \enumerate{
#'         \item{In case \code{ep$exec.procedure="Optimization"} a scalar 
#'             target value is returned.}
#'         \item{In case \code{ep$exec.procedure="Basic-Param"} a list 
#'             with transient result data is returned.}}
#'                                                                  
#' @section ToDo:
#' \itemize{
#'    \item{Capture Error Messages in foreach and model loop}
#'    \item{Fix county plots}
#'    \item{Implement statisitcs output against opt.targets}
#'    \item{Implement normed standard deviation as target value in     
#'          Global daths & icu_cases & local deaths}
#'    \item{Implement Error message in case R0county contains county id which
#'          is not selected for simulation.}
#' }
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import data.table
#' @import rlist
#' @import lhs
#' 
#' @export
CoSMic <- function(ep, sp, iol, pspace, sim.struc, op, opt) {
   
    ##    ep  <-  exec.params
    ##    sp  <-  static.params
    ##    iol = iol
    ##    pspace=pspace
    ##    sim.struc = sim.struc
    ##    op   <-  opt.params
    ##    opt  <-  rep(0.5,6) #NULL #
    
    ## Set seed for random numbers (for replication) ---------------------------
    if ( is.null( sp$seed.base ) ) {

        ## Get current seconds since epoch from system ---------------
        set.seed(as.integer(format.POSIXct(Sys.time(), format = "%s")))
        
    } else {
        set.seed(sp$seed.base)
    }
    
    ## #########################################################################
    ## Preparations                                                            #
    ## #########################################################################

    ## States that the simulated individuals can reach ----------------
    healthy    <- 0 # Healthy
    inf_noncon <- 1 # Infected, not contagious
    inf_contag <- 2 # Infected, contagious
    ill_contag <- 3 # Ill, contagious
    ill_ICU    <- 4 # Ill, ICU
    immune     <- 5 # Immune
    dead       <- 6 # Dead
    
    ## NA -----------------------------------------------------------
    missing <- -99

    ## Prepare transition probabilities information ############################
    if (sp$control_age_sex=="NONE") {
        ch_age <- c("total")
        ch_sex <- c("total")
    }
    
    if (sp$control_age_sex=="age") {
        ch_age <- c(paste(seq(0,90,5)))
        ch_sex <- c("total")
    }
    if (sp$control_age_sex=="sex") {
        ch_age <- c("total")
        ch_sex <- c("m","f")
    }
    if (sp$control_age_sex=="age_sex") {
        ch_age <- c(paste(seq(0,90,5)))
        ch_sex <- c("m","f")
    }
    
    ## Generate seed_sequence of days for realistic seed #######################

    ## We add 1 here so that ill people are generated also using
    ## the information for the day after the seed (assuming that people
    ## who are detected 1 day after seed are already ill on seed day)
    seed_date   <- sp$seed_date+1
    seed_before <- seed_date-1-sp$seed_before
    seed_seq    <- seq(from=seed_before,to=seed_date,by="days")
    #seed_seq       <- as.character(seed_seq)
    
    ## Derive dates of infections for those that are inf_cont, but are not -----
    ## yet aware about it (will be registered the next two days)
    seed_inf_cont <- seed_date+sp$cont_dur
    seed_inf_cont_seq <- seq(from=seed_date+1,to=seed_inf_cont,by="days")
    #seed_inf_cont_seq <- as.character(seed_inf_cont_seq)

    ## Derive dates of infections for those that are inf_cont, but are not -----
    ## yet aware about it (will be registered the next 3-5 days)
    seed_inf_ncont_s   <- seed_date+sp$cont_dur+1
    seed_inf_ncont_e   <- seed_date+sp$inf_dur+sp$cont_dur
    seed_inf_ncont_seq <- seq(from = seed_inf_ncont_s,
                              to   = seed_inf_ncont_e, by = "days")
    #seed_inf_ncont_seq <- as.character(seed_inf_ncont_seq)

    ## In case time_n is smaller than the max value of sp$R0change, ------------
    ## reset time_n. Otherwise a subscript out of bounds error will 
    ## be raised during R0matrix generation.
    if (is.null(sp$time_n)) {
        sp$time_n <- max(unlist(sp$R0change))+1
        save.static.params(ep, sp)
    } else if (max(unlist(sp$R0change)) > sp$time_n) {
        sp$time_n <- max(unlist(sp$R0change))+1
        save.static.params(ep, sp)
    }

    ## -------------------------------------------------------------------------
    ## Initialize the data.frame carrying the different sets of model parameters
    ## to be evaluated.
    if (sp$lhc.reload) {
        lhc <- iol$lhc.data
    } else {
        lhc <- init.lhc(pspace,sp)
    }
        
    ## -------------------------------------------------------------------------
    ## If we want to pass in parameters, opt has to be a numeric matrix --------
    if (length(opt) > 0) {

        opt.dat <- t(replicate(dim(lhc)[1],as.numeric(opt)))
        
        if ( sum(dim(lhc[,op$opt.names]) != dim(opt.dat)) > 0 ) {
            stop(paste("Something 's badly wrong !!! \n",
                       "dim(lhc[,op$opt.names]) != dim(opt.dat) \n",
                       "Please check the optimization setup !!!")   )
        }
                
        lhc[,op$opt.names] <- opt.dat
    }
    
    ## Control print of parameter space data.frame -----------------------------
    if ((ep$exec.procedure == "Basic-Param") &
        (dim(lhc)[1] <= 64)) {
        print(lhc)
    }

    ## If checkpoint reload is requested check if enough dumps are available ---
    ## if not issue warning that checkpoint data are reused among executions ---
    if (ep$cp.reload) {
        n.cp <- length(dir(ep$cp.dir,pattern="checkpoint.data-proc[0-9]*.RDS$"))

        if (n.cp == 0) {
            stop(paste("No checkpoints found in",ep$cp.dir))
        }
        
        if (n.cp < dim(lhc)[1]) {
            warning("Number of restart files is smaller than dim(lhc)[1].")
            print("Number of restart files is smaller than dim(lhc)[1].")
        }
        
    }
        
    ## Init cluster for nested parallel execution in case of optimization
    ## and parallelizaton based on a POCK-Cluster
    if ( (ep$exec.procedure  == "Optimization") &
         (ep$parallel.method == "PSOCK")        ) {
       
        ## Use nettest from VirtualGL to find a free port ----------------------
        base.port <- as.integer(system("nettest -findport",intern=TRUE))
        ##print(paste(Sys.info()["nodename"],base.port))
        
        cls <- makeClusterPSOCK(sp$iter,
                                verbose=ep$omp.cluster.dbg,
                                autoStop=FALSE,
                                outfile=paste0(getwd(),"/cls-",Sys.info()["nodename"],"-",base.port,".out"),
                                port=base.port)

        ##print(showConnections()[showConnections()[,"class"]=="sockconn",])
        registerDoParallel(cls)
        dbg.out <- clusterCall(cls,setDTthreads, threads=1)
        dbg.out <- clusterCall(cls,options, dplyr.summarise.inform=FALSE)
        
    }

    ## If R0effect is passed in as directl which means it should be changed ----
    ## per R0change, determine in which kind of region it should change,    ----
    ## either states or Nuts2.                                              ----
    if ( pspace$R0effect$type %in% c("directl","distv") ) {
        ## grepl pattern with state shortcuts ----------------------------------
        np.s <- paste(iol[["states"]]$Shortcut,collapse="|")
        ## grepl pattern with nuts2 shortcuts ----------------------------------
        np.n <- paste(iol[["counties"]]$Nuts2,collapse="|")

        ## In which region R0effect changes ------------------------------------
        if (sum(grepl(np.s,names(lhc))) > 0) {
            ## If state shortcuts are found in names(lhc) ----------------
            R0effect.region <- "state"         
        } else if (sum(grepl(np.n,names(lhc))) > 0) {
            ## If Nuts2 shortcuts are found in names(lhc) ----------------
            R0effect.region <- "nuts2"
        } else {
            stop(paste("Names in the lhc determined based on pspace are:\n",
                       paste(names(lhc),collapse=","),"\n",
                       "No state or Nuts2 shortcuts could be found.\n"))
        }
        print(paste("R0effect will change within",R0effect.region))
    }
            
################################################################################
###                                                                            #  
### Simulation model                                                           # 
###                                                                            # 
################################################################################
    
    fer <- foreach (it.ss=seq(1,dim(lhc)[1],1),
                    .packages=c("dplyr","data.table")
                    ) %dopar%
        {

            #print(showConnections()[showConnections()[,"class"]=="sockconn",])
            
            ## Threadded paralellism does not help to gain speedup 
            ## on the systems tested so we deactivate it 
            setDTthreads(threads=1)
            
            ## Set Random seed to seed.base and see whether all iterations are equal
            if (sp$seed.in.inner.loop) {
                set.seed(sp$seed.base)
            } else {
                ## Set different seeds in  all iterations to make iterations
                ## executed differently by different MPI processes.
                set.seed(as.integer(
                    as.numeric(format.POSIXct(Sys.time(), format = "%s"))/it.ss))
            }
            
            ## Chance of moving into intensive care unit when ill (per day)
            target_icu_risk <- iol[["trans_pr"]]$icu_risk[iol[["trans_pr"]]$age%in%ch_age&
                                                          iol[["trans_pr"]]$sex%in%ch_sex]
            icu_risk <- 1-(1-target_icu_risk)^(1/sp$ill_dur)
            
            ## Put into data frame
            if (sp$control_age_sex=="age") {
                icu_risk <- rep(icu_risk,2)
            }
            if (sp$control_age_sex=="sex") {
                icu_risk <- c(rep(icu_risk[1],length(seq(0,90,5))),
                              rep(icu_risk[2],length(seq(0,90,5))))
            }
            
            icu_risk <- data.frame(age=iol[["trans_pr"]]$age_gr[-which(iol[["trans_pr"]]$sex=="total"|
                                                                       iol[["trans_pr"]]$age=="total")],
                                   sex=iol[["trans_pr"]]$sex[-which(iol[["trans_pr"]]$sex=="total"|
                                                                    iol[["trans_pr"]]$age=="total")],
                                   icu_risk=icu_risk)
            
            ## ICU risk per day ------------------------------------------------
            icu_risk_list <- list()
            for (ill_days in 1:sp$ill_dur) {
                tobind <- icu_risk
                tobind$dur <- ill_days
                tobind$icu_risk <- tobind$icu_risk*sp$icu_per_day[ill_days]
                icu_risk_list[[ill_days]] <- tobind
            }    
            icu_risk <- data.frame(rbindlist(icu_risk_list))
            
            ## Survival chance of sick but not in ICU (per day) ----------------
            target_surv_ill <- iol[["trans_pr"]]$surv_ill[iol[["trans_pr"]]$age%in%ch_age&
                                                          iol[["trans_pr"]]$sex%in%ch_sex]
            surv_ill <- target_surv_ill^(1/sp$ill_dur) 
            
            ## Put into data frame ---------------------------------------------
            if (sp$control_age_sex=="age") {
                surv_ill <- rep(surv_ill,2)
            }
            if (sp$control_age_sex=="sex") {
                surv_ill <- c(rep(surv_ill[1],length(seq(0,90,5))),
                              rep(surv_ill[2],length(seq(0,90,5))))
            }
            surv_ill <- data.frame(age=iol[["trans_pr"]]$age_gr[-which(iol[["trans_pr"]]$sex=="total"|
                                                                       iol[["trans_pr"]]$age=="total")],
                                   sex=iol[["trans_pr"]]$sex[-which(iol[["trans_pr"]]$sex=="total"|
                                                                    iol[["trans_pr"]]$age=="total")],
                                   surv_ill=surv_ill)
            
            ## Survival chance when at IC unit (per day) -----------------------
            target_surv_ICU <- iol[["trans_pr"]]$surv_icu[iol[["trans_pr"]]$age%in%ch_age&
                                                          iol[["trans_pr"]]$sex%in%ch_sex]
            surv_ICU <- target_surv_ICU^(1/lhc[it.ss,"icu_dur"])
            
            ## Put into data frame
            if (sp$control_age_sex=="age") {
                surv_ICU <- rep(surv_ICU,2)
            }
            if (sp$control_age_sex=="sex") {
                surv_ICU <- c(rep(surv_ICU[1],length(seq(0,90,5))),
                              rep(surv_ICU[2],length(seq(0,90,5))))
            }
            
            surv_ICU <- data.frame(age=iol[["trans_pr"]]$age_gr[-which(iol[["trans_pr"]]$sex=="total"|
                                                                       iol[["trans_pr"]]$age=="total")],
                                   sex=iol[["trans_pr"]]$sex[-which(iol[["trans_pr"]]$sex=="total"|
                                                                    iol[["trans_pr"]]$age=="total")],
                                   surv_ICU=surv_ICU)
            
            ## Random sampling -------------------------------------------------

            if(sp$sim_pop=="random") {
                
                ## Generate sample ---------------------------------------------
                rownumbers <- 1:dim(sim.struc[["pop"]])[1]
                rownumbers <- sample(rownumbers, size=lhc[it.ss,"sam_size"],
                                     prob=sim.struc[["pop"]]$total,rep=T)
                
                ## Build sample data frame -------------------------------------
                sim <- sim.struc[["pop"]][rownumbers,c("dist_id","sex","age_gr")]
                
            }
            
            ## Proportional ----------------------------------------------------
            
            if(sp$sim_pop=="proportional") {

                ## Rescale proportionally
                sim.struc[["pop"]]$total <- round(sim.struc[["pop"]]$total/
                                                  sum(sim.struc[["pop"]]$total) *
                                                  lhc[it.ss,"sam_size"])
                
                ## Build sample data frame
                sim <- as.data.frame(lapply(sim.struc[["pop"]], rep, sim.struc[["pop"]]$total))
                
                ## Restrict variables
                sim <- sim[,c("dist_id","sex","age_gr")]
                
            }
    
            ## Expand to cover all time steps -----------------------------
            ## Ralf: This is only possible for sam_size*time_n < 2^30.
            ##       I.e. with time_n = 100 for sam_size_n < 20 Million.
            ##       R is not able to address dataframes with more entries.
            ## sim[,paste0("t",1:time_n)] <- missing
            
            ## Expand to cover two time steps and do a Flip-Flop
            sim[,paste0("t",1:2)] <- missing
    
            ## States for first time period
            sim$t1 <- healthy
    
            ## Duration variable
            sim$d <- 1

### Seed infections ###########################################################

            ## Define initially infected population (random) -------------------
    
            if(sp$seed_infections=="random") {
    
                ## Reshuffle cases
                sim <- sim[sample(1:dim(sim)[1]),]
                
                ## Assign infection
                sim[1:sp$ini_infected,"t1"] <- inf_contag
                
            }
            
            ## Define initially infected population (random county)
            
            if(sp$seed_infections=="county") {
                
                ## Sample random county
                rcounty <- sample(unique(sim$dist_id),1)
                
                ## Number of individuals in county
                ncounty <- sum(sim$dist_id==rcounty)
      
                ## Infection probability for county
                pcounty <- min(sp$ini_infected/ncounty,1)
      
                ## Sample infections
                iol[["seed"]] <- sample(c(healthy,inf_contag),
                                        size=ncounty,
                                        rep=T,
                                        prob=c(1-pcounty,pcounty))
      
                ## Assign infections
                sim[sim$dist_id==rcounty,"t1"] <- iol[["seed"]]
                
            }
    
            ## Realistic seed --------------------------------------------------
    
            if(sp$seed_infections=="data") {
      
                ## Restrict counties
                ## If there are no counties in the model, the first dist_id is used
                if(all(sim.struc$counties=="NONE")) iol[["seed"]]$distid <- unique(sim$dist_id)[1]
                if (all(sim.struc$counties=="NONE")) iol[["seed_dea"]]$distid <- unique(sim$dist_id)[1]
      
                ## Derive ids of all counties in the model
                sim_counties <- unique(sim$dist_id)
                ## Select data for these districts
                iol[["seed"]] <- iol[["seed"]][iol[["seed"]]$distid%in%sim_counties,]
                iol[["seed_dea"]] <- iol[["seed_dea"]][iol[["seed_dea"]]$distid%in%sim_counties,]
                
                ## If there are no counties in the model, data for whole Germany is aggregated
                if(all(sim.struc$counties=="NONE")) iol[["seed"]] <-
                                              aggregate(cases~distid+
                                                            date,data=iol[["seed"]],sum)
                
                if(all(sim.struc$counties=="NONE")) iol[["seed_dea"]] <-
                                              aggregate(deaths~distid+
                                                            date,data=iol[["seed_dea"]],sum)
                
                ## Define those that are in the states inf_noncon, inf_contag,
                ## and ill_contag
                ## Select dates prior to start date to derive ill persons
                seed_ill <- iol[["seed"]][iol[["seed"]]$date%in%seed_seq,]
                ## Select dates after start date to derive infected, contagious persons
                seed_inf_cont <- iol[["seed"]][iol[["seed"]]$date%in%seed_inf_cont_seq,] 
                ## Select dates after start date to derive infected, not contagious persons
                seed_inf_ncont <- iol[["seed"]][iol[["seed"]]$date%in%seed_inf_ncont_seq,]
                ## Define those who are dead
                ## Select dates up until start date
                ## Derive dates for deaths (up until start date)
                ## Need to move start data back by one, as we moved it forward one
                ## day further above
                seed_d_seq <- seq(from = as.Date(min(iol[["seed_dea"]]$date)),
                                  to   = seed_date-1, by = "days")
                seed_d_seq<- as.character(seed_d_seq)
                
                ## Select dataset up until start date
                iol[["seed_dea"]] <- iol[["seed_dea"]][which(iol[["seed_dea"]]$date%in%seed_d_seq),]
                
                ## Derive durations of illness from cases recorded prior to the start date
                ## (based on the assumption that in order to be detected as a case, persons
                ## need to be at least ill for a day)
                ## Ill
                seed_ill$dur <- as.numeric(seed_date - as.Date(paste(seed_ill$date)))+1
                
                ## Exclude all cases that were recorded before the illness duration period, 
                ## and all cases below 1
                seed_ill <- seed_ill[which(seed_ill$dur<=sp$ill_dur&seed_ill$cases>0),]
                
                ## Derive durations for infected, contagious
                seed_inf_cont$dur <- as.numeric(seed_date+sp$cont_dur+1-
                                                as.Date(paste(seed_inf_cont$date)))
                
                ## Exclude all cases below 1
                seed_inf_cont <- seed_inf_cont[which(seed_inf_cont$cases>0),]
                
                ## Derive durations for infected, not contagious
                seed_inf_ncont$dur <- as.numeric(seed_date +
                                                 sp$inf_dur+sp$cont_dur+1-
                                                 as.Date(paste(seed_inf_ncont$date)))
                ## Exclude all cases below 1
                seed_inf_ncont <- seed_inf_ncont[which(seed_inf_ncont$cases>0),]

                ## Aggregate ill
                seed_ill <- aggregate(cases~distid+dur,data=seed_ill,sum)
                
                ## Aggregate infected, contagious
                seed_inf_cont <- aggregate(cases~distid+dur,data=seed_inf_cont,sum)
                
                ## Aggregate infected, not contagious
                seed_inf_ncont <- aggregate(cases~distid+dur,data=seed_inf_ncont,sum)
                
                ## Aggregate deaths
                seed_dth <- aggregate(deaths~distid,data=iol[["seed_dea"]],sum)
                
                ## Proportional to sample size, ill ---------------------------------------
                ##seed_ill$cases <- ceiling(seed_ill$cases*sam_prop)
                seed_ill$cases <- apply(seed_ill,1,
                                        function(x){
                                            ceiling(x[3] *  sp$sam_prop.ps[as.integer(x[1]/1000)])})
                
                ## Proportional to sample size, infected, contagious ----------------------
                ## seed_inf_cont$cases <- ceiling(seed_inf_cont$cases*sam_prop)
                seed_inf_cont$cases <- apply(seed_inf_cont,1,
                                             function(x){
                                                 ceiling(x[3] * sp$sam_prop.ps[as.integer(x[1]/1000)])})
                
                ## Proportional to sample size, infected, not contagious 
                ##seed_inf_ncont$cases <- ceiling(seed_inf_ncont$cases*sam_prop)
                seed_inf_ncont$cases <- apply(seed_inf_ncont,1,
                                              function(x){
                                                  ceiling(x[3] *  sp$sam_prop.ps[as.integer(x[1]/1000)])})
                
                ## Proportional to sample size, deaths
                ## seed_dth$deaths <- ceiling(seed_dth$deaths*sam_prop)
                seed_dth$deaths <- apply(seed_dth,1,
                                         function(x){
                                             ceiling(x[2] *  sp$sam_prop.ps[as.integer(x[1]/1000)])})

                ## Loop over counties ------------------------------------------
                for(county in sim_counties) {
                    
                    ## Detect row numbers of district  
                    rownumbers <- which(sim$dist_id==county)
                    
                    ## Derive durations for ill
                    se_il <- seed_ill[seed_ill$distid==county,]
                    il_d <- rep(se_il$dur,se_il$cases)
                    
                    ## Derive durations for infected, contagious 
                    se_inf_c <- seed_inf_cont[seed_inf_cont$distid==county,]
                    inf_c_d <- rep(se_inf_c$dur,se_inf_c$cases)
                    
                    ## Derive durations for infected, not contagious 
                    se_inf_nc <- seed_inf_ncont[seed_inf_ncont$distid==county,]
                    inf_nc_d <- rep(se_inf_nc$dur,se_inf_nc$cases)
              
                    ## Derive total ill/infected, contagious/infected, not contagious
                    inf_ill <- max(sum(seed_ill[seed_ill$distid==county,"cases"]),0)
                    inf_cont <- max(sum(seed_inf_cont[seed_inf_cont$distid==county,
                                                      "cases"]),0)
                    inf_ncont <- max(sum(seed_inf_ncont[seed_inf_ncont$distid==county,
                                                        "cases"]),0)
                    
                    ## Derive total nuber of dead
                    inf_dth <- max(sum(seed_dth[seed_dth$distid==county,
                                                "deaths"]),0)

                    if(length(rownumbers)<c(inf_ill+inf_cont+inf_ncont+inf_dth)) 
                        stop("Number of infected and dead is larger than population size")
                    
                    ## Determine rownumbers of ill persons
                    rownumbers_ill <- 0
                    rownumbers_left1 <- rownumbers
                    if (inf_ill >0) {
                        rownumbers_ill <- sample(rownumbers,size=inf_ill)         
                        rownumbers_left1 <- rownumbers[-which(rownumbers%in%rownumbers_ill)]
                    }
                    
                    ## Determine rownumbers of infected, contagious persons
                    rownumbers_cont <- 0
                    rownumbers_left2 <- rownumbers_left1
                    if (inf_cont >0) {
                        rownumbers_cont <- sample(rownumbers_left1,size=inf_cont)
                        rownumbers_left2 <- rownumbers_left1[-which(rownumbers_left1%in%rownumbers_cont)]
                    }
                    ## Determine rownumbers of infected, contagious persons
                    rownumbers_ncont <- 0
                    rownumbers_left3 <- rownumbers_left2
                    if (inf_ncont >0) {
                        rownumbers_ncont <- sample(rownumbers_left2,size=inf_ncont)
                        rownumbers_left3 <- rownumbers_left2[-which(rownumbers_left2%in%rownumbers_ncont)]
                    }
                    ## Determine rownumbers for death persons
                    rownumbers_dea <- 0
                    if (inf_dth>0) {
                        rownumbers_dea <- sample(rownumbers_left3,size=inf_dth)
                    }
                    ## Assign states and durations
                    sim[rownumbers_ill,"t1"] <- ill_contag
                    sim[rownumbers_ill,"d"] <- il_d
                    sim[rownumbers_cont,"t1"] <- inf_contag
                    sim[rownumbers_cont,"d"] <- inf_c_d
                    sim[rownumbers_ncont,"t1"] <- inf_noncon
                    sim[rownumbers_ncont,"d"] <- inf_nc_d
                    sim[rownumbers_dea,"t1"] <- dead
                }
            } 
            ## Set county vector to "fake" county if no county structure
            if(all(sim.struc$counties=="NONE")) sim.struc$counties <- unique(sim$dist_id)[1]
            
### Generate R0 matrix: Daily R0 ##############################################
            
            if (sp$import_R0_matrix==TRUE) {
                R0_raw1 <- iol[["R0_raw"]][match(sim.struc$counties,iol[["R0_raw"]]$dist_id),]
                sp$R0change <- as.matrix(R0_raw1[,-1])
            }
            
            ## Daily R0
            ## R0_daily <- sp$R0_force*(R0/(sp$cont_dur+sp$ill_dur*sp$less_contagious)) +
            ##   (1-sp$R0_force)*R0/(sp$cont_dur+sp$ill_dur)
            R0_daily <- sp$R0_force*(lhc[it.ss,"R0"]/(sp$cont_dur+sp$ill_dur*sp$less_contagious)) +
                (1-sp$R0_force)*lhc[it.ss,"R0"]/(sp$cont_dur+sp$ill_dur)
            
            ## Matrix filled with R0_daily
            R0matrix <- matrix(data=R0_daily,nrow=length(sim.struc$counties),ncol=sp$time_n-1)
            
            ## Rownames = county names
            rownames(R0matrix) <- paste(sim.struc$counties)
            
            ## Get number of changes
            if(all(sp$R0change!="NONE") & 
               !is.matrix(sp$R0change)) nchange <- length(sp$R0change) else nchange <- 0
            
            ## Edit R0 values if sp$R0change = list
            if(is.matrix(sp$R0change)) R0matrix <- sp$R0change * R0matrix
            
            ## Edit R0 values
            if(nchange>0) {
                
                ## Loop over changes
                for(change in 1:nchange) {
                    
                    ## Get all affected counties
                    getcounties <- sp$R0county[[change]]
                    if(all(getcounties=="ALL")) getcounties <- paste(sim.struc$counties)

                    ## Get relative change ----------------------------------------

                    if ( is.null(pspace$R0effect) ) {
                        stop("R0effect is missing in pspace!")
                    }

                    ## If we have a directl type for R0effect ------------------
                    if ( pspace$R0effect$type %in% c("directl","distv") ) {

                        ## If R0effect changes per state per change ------------
                        if (R0effect.region == "state") {
                            
                            getchange <- as.numeric(t(lhc[it.ss,
                                                          paste0(iol[["states"]][
                                                              as.character(as.integer(as.integer(
                                                                  getcounties)/1000)),"Shortcut"],change)]))
                        }

                        ## If R0effect changes per Nuts2 region per change -----
                        if ( R0effect.region == "nuts2" ) {

                            getchange <- as.numeric(t(lhc[it.ss,
                                                          paste0(iol[["counties"]][
                                                              iol[["counties"]]$dist_id%in%getcounties,
                                                              "Nuts2"],change)]))
                        } 
                       
                    } else if ( pspace$R0effect$type == "directv" ) {
                        ## Per change ------------------------------------------
                        getchange <- lhc[it.ss,paste0("R0effect",change)]
                        
                    } else {
                        ## Globally constant -----------------------------------
                        getchange <- lhc[it.ss,"R0effect"]
                    }
                    
                    ## Get affected time steps
                    gettime <- sp$R0change[[change]]
                    gettime <- seq(gettime[1],gettime[2],by=1)

                    ## Change R0 values
                    R0matrix[getcounties,gettime] <- R0matrix[getcounties,gettime]*getchange
                    
                }
            }
            
            ## Delay and smooth
            if(sp$R0delay) {
                
                ## Apply function row-wise
                R0matrix <- t(apply(R0matrix,1,attenuate,
                                    steps=sp$R0delay_days,type=sp$R0delay_type))
            }
            
### Objects for additional results ############################################
### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            ## Object for results on infections
            inf_cases <- matrix(data=NA,ncol=sp$time_n,nrow=length(sim.struc$counties))
            rownames(inf_cases) <- paste(sim.struc$counties)
            
            ## Objects for ICU cases
            icu_cases <- inf_cases

            ## Objects to store total daily numbers ----------------------
            healthy_cases    <- inf_cases
            inf_noncon_cases <- inf_cases
            inf_contag_cases <- inf_cases
            ill_contag_cases <- inf_cases
            ill_ICU_cases    <- inf_cases   
            immune_cases     <- inf_cases    
            dead_cases       <- inf_cases      
            dead_cases_bICU  <- inf_cases      
            mod_inf_cases    <- inf_cases
            org_noncon_cases <- inf_cases
            
            ## Object to track lockdowns
            lockdowns <- matrix(data=NA,ncol=sp$time_n,nrow=1)
            
            ## Get starting values for infections
            start_values_tot <- by(sim$t1,sim$dist_id,
                                   function(x) sum(table(x)[c(as.character(inf_contag),
                                                              as.character(inf_noncon),
                                                              as.character(ill_contag))]))
            start_values_tot[is.na(start_values_tot)] <- 0
            
            ## Assign starting values ----------------
            inf_cases[,1]        <- c(start_values_tot)
            icu_cases[,1]        <- 0
            
            start_values <- by(sim$t1,sim$dist_id,
                               function(x) table(x)[paste(healthy)])
            start_values[is.na(start_values)] <- 0
            healthy_cases[,1]    <- c(start_values)
            
            start_values <- by(sim$t1,sim$dist_id,
                               function(x) table(x)[paste(inf_noncon)])
            start_values[is.na(start_values)] <- 0
            inf_noncon_cases[,1]    <- c(start_values)
            
            start_values <- by(sim$t1,sim$dist_id,
                               function(x) table(x)[paste(inf_contag)])
            start_values[is.na(start_values)] <- 0
            inf_contag_cases[,1] <- c(start_values)
            
            start_values <- by(sim$t1,sim$dist_id, 
                               function(x) table(x)[paste(ill_contag)])
            start_values[is.na(start_values)] <- 0
            ill_contag_cases[,1]    <- c(start_values)
            
            start_values <- by(sim$t1,sim$dist_id,
                               function(x) table(x)[paste(dead)])
            start_values[is.na(start_values)] <- 0
            dead_cases[,1]    <- c(start_values)
            ill_ICU_cases[,1]    <- 0   
            immune_cases[,1]     <- 0    
            dead_cases_bICU[,1]  <- 0    
            mod_inf_cases[,1]    <- 0
            org_noncon_cases[,1] <- 0
            
            max.case.date <- max(as.Date(iol[["seed"]]$date))

            ## If Checkpoint reload is requested ----------------------------
            if (ep$cp.reload) {
                cp.data <- readRDS(paste0(ep$cp.dir,"/checkpoint.data-proc",
                                          ((it.ss-1)%%n.cp)+1,
                                          ".RDS"))

                time_n <- dim(cp.data$healthy_cases)[2]
                
                healthy_cases[,c(1:time_n)]      = cp.data$healthy_cases   
                inf_noncon_cases[,c(1:time_n)]   = cp.data$inf_noncon_cases
                inf_contag_cases[,c(1:time_n)]   = cp.data$inf_contag_cases
                ill_contag_cases[,c(1:time_n)]   = cp.data$ill_contag_cases
                ill_ICU_cases[,c(1:time_n)]      = cp.data$ill_ICU_cases   
                immune_cases[,c(1:time_n)]       = cp.data$immune_cases    
                dead_cases[,c(1:time_n)]         = cp.data$dead_cases      
                inf_cases[,c(1:time_n)]          = cp.data$inf_cases       
                icu_cases[,c(1:time_n)]          = cp.data$icu_cases       
                dead_cases_bICU[,c(1:time_n)]    = cp.data$dead_cases_bICU 
                mod_inf_cases[,c(1:time_n)]      = cp.data$mod_inf_cases   
                org_noncon_cases[,c(1:time_n)]   = cp.data$org_noncon_cases
                lockdowns[,c(1:time_n)]          = cp.data$lockdowns
                
                sim = cp.data$sim

                start.time <- time_n + 1

                print(paste("Checkpoint reloaded. start.time of simulation",
                            "loop is set to",start.time))
                      
            } else {
                start.time <- 2
            }
            
            
################################################################################
### Simulation Loop ############################################################
            print("Simulation Loop")
            
            for(timestep in start.time:sp$time_n) {

                print(timestep)
                
                ## Cat for each day
                ## cat(".")
    
                ## Get previous state variable
                prev_state <- paste0("t",timestep-1)
                
                ## Get state to generate
                up_state <- paste0("t",timestep)

                ## Ralf Schneider: To save one %in% Operation the logic is changed to:
                ## First reduce sim to tmp and subset and then check whether all are immune
                ## or dead
                
                ## Check if anything to do (all dead/immune?)
                ##if(all( sim[,prev_state] %in% c(immune,dead) )) {
                ##  
                ##  ## If all dead/immune: Keep everything unchanged
                ##  sim[,up_state] <- sim[,prev_state]
                ##  cumulative[timestep] <- cumulative[timestep-1] 
                ##  next  
                ##}
                ##
                ## Select variables of interest 
                ##tmp <- sim[,c("dist_id","age_gr","d",prev_state,up_state)]
                ##
                ## Subset of the living and NOT immune 
                ##tmp <- tmp[! tmp[,prev_state] %in% c(immune,dead),]

                ## Select variables of interest 
                tmp <- sim[,c("dist_id","age_gr","sex","d",prev_state,up_state)] 
                
                ## Subset of the living and NOT immune
                ##tmp.dplyr <- tmp %>% filter(tmp[,prev_state] != immune & tmp[,prev_state] != dead)
                ##print(dim(tmp.dplyr))
                
                ##tmp <- tmp[! tmp[,prev_state] %in% c(immune,dead),]
                ##print(dim(tmp))
                tmp <- tmp %>% filter(tmp[,prev_state] != immune & tmp[,prev_state] != dead)
                
                ## Check if anything to do (all dead/immune?)
                ##if(all( sim[,prev_state] %in% c(immune,dead) )) {
                if (dim(tmp)[1] == 0) {
                    
                    ## If all dead/immune: Keep everything unchanged
                    sim[,up_state] <- sim[,prev_state]
                    cumulative[timestep] <- cumulative[timestep-1] 
                    next  
                }
                
                ## New value of counting variable
                tmp$d_new <- missing

                ################################################################
                ## State: Healthy, susceptible #################################
                
                ## Susceptible population by county ============================

                ## Who is healthy ---------------------------
                tmp$count   <- tmp[,prev_state] == healthy

                ## How many susceptibles per state ----------
                susceptible <- group_by(tmp, dist_id) %>% summarise(count = sum(count))
      
                ## Contagious population by county ==============================

                ## Who is contagious ------------------------
                part1      <-   tmp[,prev_state] ==  inf_contag
                part2      <- ( tmp[,prev_state] ==  ill_contag ) * sp$less_contagious
                tmp$count  <- part1 + part2
                
                ## How many contagious per state ------------
                contagious <- group_by(tmp, dist_id) %>% summarise(count = sum(count))
            
                ## Check if susceptible and contagious are around ===============
                at_risk      <- sum(susceptible$count)
                n_contagious <- sum(contagious$count)
                
                ## If infections (potentially) possible =========================
                if(at_risk>0 & n_contagious>0) {
                    
                    ## Expected infections by county ----------
                    exp_infect <- contagious$count * R0matrix[paste(contagious$dist_id),
                                                              timestep-1]
        
                    ## If endogenous lockdown ----------------------------------
                    if(sp$endogenous_lockdown) {
                        
                        ## Indices of days to be checked ----
                        checkdays <- (timestep-1):(timestep-sp$lockdown_days)
                        
                        ## Remove invalid indices
                        checkdays <- checkdays[!checkdays<1]
                        
                        ## Get daily counts of new cases
                        checkdays <- colSums(inf_cases)[checkdays]
                        
                        ## Check if any are equal/above threshold
                        checkdays <- any(checkdays>=sp$lockdown_threshold)
                        
                        ## If any daily counts above threshold: Lockdown
                        if(checkdays) {
                            
                            ## Apply lockdown effect
                            exp_infect <- exp_infect * sp$lockdown_effect
                            
                            ## Keep track of lockdowns
                            lockdowns[timestep] <- 1
                        } else lockdowns[timestep] <- 0
                        
                    } else lockdowns[timestep] <- 0
                    
                    ## Rescale connectivity (if not all counties have contagious) 
                    connect <- iol[["connect_work"]][paste0(contagious$dist_id),
                                                     paste0(contagious$dist_id)]
                    connect <- as.matrix(connect)
                    connect <- apply(connect,1,function(x) x/sum(x))
                    connect   <- t(connect)
                    
                    ## Spread expected infections by connectivity
                    between_weight <- 1 - lhc[it.ss,"w_int"]
                    
                    ## If lockdown in place, effect ---------
                    if(sp$endogenous_lockdown) between_weight <- between_weight * sp$lockdown_connect
                    
                    within_weight <- 1 - between_weight
                    
                    ## Apply spatial weights ----------------
                    exp_infect <- within_weight  * exp_infect +
                        between_weight * c(exp_infect%*%connect)
                    
                    ## Risk of infection by county: Denominator 
                    if(!sp$immune_stop) {
                        
                        ## Only take susceptible if immune people do not reduce spread
                        denominator <- susceptible$count
                        
                    } else {
                        
                        ## Select variables of interest -----
                        denominator <- setDT(sim[,c("dist_id",prev_state)])
                        
                        ## Select counties with susceptible -
                        denominator <- denominator %>% filter(dist_id %in% susceptible$dist_id)
                        
                        ## Base denominator on all living ---
                        denominator$count <-
                            denominator[,get(prev_state)]%in%c(healthy,immune,
                                                               inf_noncon,inf_contag,
                                                               ill_contag)
                        
                        ## Aggregate to county level --------
                        denominator <- group_by(denominator, dist_id) %>% summarise(count = sum(count))
                        denominator <- denominator$count
                    }
                    
                    ## Risk of infection by county: Calculate
                    risk <- exp_infect/denominator
                    
                    ## Capture some exceptions and avoid weird results
                    risk[is.infinite(risk)] <- 0
                    risk[is.nan(risk)] <- 0
                    risk[risk>1] <- 1
                    
                    ## Weight risk (because the healthy move)
                    risk <- lhc[it.ss,"w_int"] * risk +
                        (1-lhc[it.ss,"w_int"]) * c(risk%*%connect)
                    
                    ## Probablities for each individual at risk (by county)
                    probs <- match(tmp[tmp[,prev_state]==0,"dist_id"],contagious$dist_id)
                    probs <- risk[probs]
                    
                    ## Who gets sick? Randomized ------------
                    sick <- as.numeric(runif(at_risk)<=probs)
                    initial.sick  <- sum(sick)
                    
                    ## Count new infections by county -------
                    cbtmp<-data.table(sick,dist_id=tmp[tmp[,prev_state]==0,"dist_id"])
                    
                    case_count<-group_by(cbtmp,dist_id) %>% summarize(count = sum(sick))
                    
                    case_count[is.na(case_count)] <- 0
                    
                    ## In case not all counties have cases spread them out -----
                    ## to full list                                        -----
                    final_count <- numeric(length(sim.struc$counties))
                    
                    final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$count
                    
                    inf_cases[,timestep] <- final_count
                    
                    org_noncon_cases[,timestep] <- final_count
                    
                    ## In case we ==============================================
                    ## - have set w.obs > 0. and                             ===
                    ## - are in the date range in which reported case        ===
                    ##   numbers exist and                                   ===
                    ## - people have gotten sick,                            ===
                    ## the infections are redistributed accross the counties ===
                    ## according to w.obs by county                          ===

                    ## Since infections are reported with 7 days delay we ---
                    ## look 7 days ahead                                  ---
                    target.date <- seed_date + timestep-1 + 7
                    
                    if ( (lhc[it.ss,"w.obs"] > 0.) &
                         (target.date < max.case.date) &
                         (initial.sick > 0)                                      ) {
                        
                        ## Get the number of infections per county from seed  --
                        target.inf  <- iol[["seed"]][
                            as.Date(iol[["seed"]]$date) == target.date ,]
                        
                        ## In case not all counties have observed cases   ------
                        ## spread them out to full list                   ------
                        final_count <- numeric(length(sim.struc$counties))
                        
                        final_count[match(target.inf$distid,sim.struc$counties)] <- target.inf$cases
                        
                        target.inf <- final_count
                        names(target.inf) <-  rownames(inf_cases)
                        
                        ## Calculate target proportion per county in case ---
                        ## it is something to distribute                  ---
                        if ( sum(target.inf) > 0 ) {
                            
                            if ( lhc[it.ss,"w.obs.by.state"] == 1 ) {
                                
                                state.ids <- unique(as.integer(sim.struc$counties/1000))
                                
                                mod.inf <- numeric(length(sim.struc$counties))
                                names(mod.inf) <- rownames(inf_cases)
                                
                                for ( st in state.ids ) {
                                    
                                    ## Proportion of cases per county ----------
                                    select.counties <- which(as.integer(sim.struc$counties/1000)==st)
                                    
                                    prop.inf_cases  <- inf_cases[select.counties,timestep]
                                    initial.sick.per.state <-  sum(prop.inf_cases)
                                    
                                    if ( initial.sick.per.state > 0 ) { 
                                        prop.inf_cases <-  prop.inf_cases / initial.sick.per.state 
                                    } else {
                                        prop.inf_cases <- 0
                                    }
                                    
                                    prop.target.inf <- target.inf[select.counties]
                                    target.sick.per.state <- sum(prop.target.inf)
                                    
                                    if ( target.sick.per.state > 0 ) {
                                        prop.target.inf <- prop.target.inf / target.sick.per.state
                                    } else {
                                        prop.target.inf <- 0
                                    }
                                    
                                    prop.target <- lhc[it.ss,"w.obs"] * prop.target.inf +
                                        (1- lhc[it.ss,"w.obs"]) * prop.inf_cases
                                    
                                    ## Calculate how many cases to shift per county ---
                                    mod.inf[select.counties] <-
                                        round(prop.target*initial.sick.per.state) -
                                        inf_cases[select.counties,timestep]
                                }
                                
                            } else {
                                
                                ## Proportion of cases per county -------
                                prop.inf_cases  <- inf_cases[,timestep] / initial.sick
                                
                                prop.target <- lhc[it.ss,"w.obs"] * ( target.inf / sum(target.inf) ) +
                                    (1- lhc[it.ss,"w.obs"]) * prop.inf_cases
                                
                                ## Calculate how many cases to shift per county ---
                                mod.inf <- round(prop.target*initial.sick)-inf_cases[,timestep]
                            }
                            
                            ## Create data.table with columns "dist.id" and "sick" ------
                            sick.dt <- data.table(dist.id=tmp[tmp[,prev_state]==0,"dist_id"],sick=sick)
                            
                            ## Shift cases according to mod.inf -------------------------
                            sick.new <- sick.dt %>% group_by(dist.id) %>%
                                group_modify(
                                    function(x,y,mod.i){
                                        
                                        mi <- mod.i[as.character(y)]
                                        
                                        if ( mi < 0 ) {
                                            to.modify <- which(x==1)
                                            if ( length(to.modify) > 0) {
                                                if (length(to.modify) > abs(mi)) {
                                                    x[sample(to.modify,abs(mi)),] <- 0
                                                } else {
                                                    x[to.modify,] <- 0
                                                }
                                            }
                                        }
                                        
                                        if ( mi > 0 ) {
                                            to.modify <- which(x==0)
                                            if (length(to.modify) > 0) {
                                                if (length(to.modify) > mi) {
                                                    x[sample(to.modify, mi ),] <- 1
                                                } else {
                                                    x[to.modify,] <- 1
                                                }
                                            }
                                        }
                                        
                                        return(x)
                                    },mod.inf)
                            
                            ## Count new infections by county ------------------
                            case_count<-group_by(sick.new,dist.id) %>% summarize(count = sum(sick))
                            
                            case_count[is.na(case_count)] <- 0
                            
                            ## In case not all counties have cases spread ------
                            ## them out to full list                      ------
                            final_count <- numeric(length(sim.struc$counties))
                            
                            final_count[match(case_count$dist.id,sim.struc$counties)] <- case_count$count
                            
                            inf_cases[,timestep] <- final_count
                            
                            mod_inf_cases[,timestep] <- mod.inf
                            
                        } else {
                            ## DEBUG printout >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                            print("No cases shifted since sum(target.inf) == 0")
                            ## DEBUG printout <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        }
                        
                    } else {
                        
                        mod_inf_cases[,timestep] <- 0
                    }
                    
                    ## Recode to match state codes ----------
                    sick[sick == 0]  <- healthy
                    sick[sick == 1]  <- inf_noncon
                    
                    ## Assign state -------------------------
                    tmp[tmp[,prev_state]==healthy,up_state] <- sick
                    
                    ## Update days (If moved to sick: First day in sickness)
                    tmp[tmp[,prev_state]==healthy & tmp[,up_state]==inf_noncon,"d_new"] <- 1
                    
                    ## Update days (If stayed healthy: Current count plus 1)
                    tmp[tmp[,prev_state]==healthy & tmp[,up_state]==healthy,"d_new"] <- 
                        tmp[tmp[,prev_state]==healthy & tmp[,up_state]==healthy,"d"] + 1
                    
                }
                
                ## If susceptible are left but no contagious individuals ========
                if(at_risk>0 & n_contagious==0) {
                    
                    ## Healthy stays healthy
                    tmp[tmp[,prev_state]==healthy,up_state] <- healthy
                    
                    ## Change days (also for people staying healthy)
                    tmp[tmp[,prev_state]==healthy & tmp[,up_state]==healthy,"d_new"] <- 
                        tmp[tmp[,prev_state]==healthy & tmp[,up_state]==healthy,"d"] + 1
                    
                    ## No new infections
                    inf_cases[,timestep] <- 0
                    
                }
                
                ## If no susceptible are left
                if(at_risk==0) {
                    
                    ## No new infections
                    inf_cases[,timestep] <- 0
                    
                }

                ## #############################################################
                ## State: Infected, non-contagious #############################
                
                ## If day limit not reach: Stays the same
                tmp[tmp[,prev_state]==inf_noncon & tmp[,"d"] <sp$inf_dur,up_state] <- inf_noncon
                
                ## Update day (Current count plus 1)
                tmp[tmp[,prev_state]==inf_noncon & tmp[,up_state]==inf_noncon,"d_new"] <- 
                    tmp[tmp[,prev_state]==inf_noncon & tmp[,up_state]==inf_noncon,"d"] + 1
                
                ## If day limit reached: Move to contagious
                tmp[tmp[,prev_state]==inf_noncon & tmp[,"d"] >=sp$inf_dur,up_state] <- inf_contag
                
                ## Update days (Reset counter to 1 if moved)
                tmp[tmp[,prev_state]==inf_noncon & tmp[,up_state]==inf_contag,"d_new"] <- 1

                ## #############################################################
                ## State: Infected, contagious  ################################
                
                ## If day limit not reach: Stays the same
                tmp[tmp[,prev_state]==inf_contag & tmp[,"d"] <sp$cont_dur,up_state] <- inf_contag
                
                ## Update days (current count plus 1 if state stays the same)
                tmp[tmp[,prev_state]==inf_contag & tmp[,up_state]==inf_contag,"d_new"] <- 
                    tmp[tmp[,prev_state]==inf_contag & tmp[,up_state]==inf_contag,"d"] + 1
                
                ## If day limit reached: Move to contagious
                tmp[tmp[,prev_state]==inf_contag & tmp[,"d"] >=sp$cont_dur,up_state] <- ill_contag
                
                ## Update days (reset to 1 if becoming ill/contagious)
                tmp[tmp[,prev_state]==inf_contag & tmp[,up_state]==ill_contag,"d_new"] <- 1

                ## #############################################################
                ## State: Ill, contagious
                
                ## Dying: Population at risk
                tmp$count   <- tmp[,prev_state]%in%c(ill_contag)

                ##at_risk <- aggregate(count~dist_id,data=tmp,sum)
                at_risk <- group_by(tmp, dist_id) %>% summarise(count = sum(count))
                
                at_risk <- sum(at_risk$count)
                
                ## Check if anyone can die
                if(at_risk>0) {
                    
                    ## Dying: Who dies?
                    age_sex <- tmp[tmp[,prev_state]%in%c(ill_contag),c("age_gr","sex")]
                    age_sex <- paste(age_sex$age_gr,age_sex$sex)
                    age_sex <- match(age_sex,paste(surv_ill$age,surv_ill$sex))
                    surv_ill_i <- surv_ill$surv_ill[age_sex]
                    die <- as.numeric(runif(at_risk)>=surv_ill_i)

                    ## Recode  
                    die[die == 0]  <- ill_contag
                    die[die == 1]  <- dead
                    ##die <- recode(die,'0' = ill_contag, '1' = dead)
                    
                    ## Dying: Assign state
                    tmp[tmp[,prev_state]==ill_contag,up_state] <- die
                }
                
                ## Save dead before moving to ICU ------------------------------
                cbtmp <- data.table(tmp[tmp[,up_state]==dead,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                dead_cases_bICU[,timestep] <- final_count
                
                ## Moving to ICU: Population at risk
                tmp$count <- tmp[,prev_state]%in%c(ill_contag) & !tmp[,up_state]%in%c(dead)
                
                ##at_risk <- aggregate(count~dist_id,data=tmp,sum)
                at_risk <- group_by(tmp, dist_id) %>% summarise(count = sum(count))

                at_risk <- sum(at_risk$count)
                
                ## Check if anyone can move to ICU
                if(at_risk>0) {
                    
                    ## Moving to ICU: Who moves? Get their ages, sex and duration
                    age_sex_dur <- tmp[tmp[,prev_state]%in%c(ill_contag) & !tmp[,up_state]%in%c(dead),
                                       c("age_gr","sex","d")]
                    age_sex_dur <- paste(age_sex_dur$age_gr,age_sex_dur$sex,age_sex_dur$d)
                    
                    ## Match with ages for ICU risk
                    age_sex_dur <- match(age_sex_dur,paste(icu_risk$age,icu_risk$sex,icu_risk$dur))
                    
                    ## Assign individual risks
                    icu_risk_i <- icu_risk$icu_risk[age_sex_dur]
                    
                    ## Draw random numbers
                    ICU <- as.numeric(runif(at_risk)<=icu_risk_i)
                    
                    ## Count new icu cases
                    ## case_count <- by(ICU,
                    ##                 tmp[tmp[,prev_state]==ill_contag & !tmp[,up_state]%in%c(dead),
                    ##                    "dist_id"],
                    ##                 function(x) table(x)["1"])
                    
                    cbtmp <- data.table(ICU,dist_id=tmp[tmp[,prev_state]==ill_contag &
                                                        !tmp[,up_state]%in%c(dead),"dist_id"])
                    
                    case_count <- group_by(cbtmp,dist_id) %>% summarize(count = sum(ICU))
                                        
                    case_count[is.na(case_count)] <- 0
                    ##case_count[is.na(case_count$count),count] <- 0

                    final_count <- numeric(length(sim.struc$counties))
                    
                    ##final_count[match(names(case_count),counties)] <- case_count
                    final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$count
                    
                    icu_cases[,timestep] <- final_count
                    
                    ## Recode --------------------
                    ##ICU <- recode(ICU,'0' = ill_contag, '1' = ill_ICU)
                    ICU[ICU == 0]  <- ill_contag
                    ICU[ICU == 1]  <- ill_ICU
                    
                    ## Moving to ICU: Assign state
                    tmp[tmp[,prev_state]==ill_contag & !tmp[,up_state]%in%c(dead),up_state] <- ICU
                    
                    ## Moving to ICU: Reset days 
                    tmp[tmp[,prev_state]==ill_contag & tmp[,up_state]==ill_ICU,"d_new"] <- 1
                    
                }
                
                ## Do not increase count if no one moves
                if(at_risk==0) {
                    icu_cases[,timestep] <- 0
                }
                
                ## Staying ill/contagious (if not dies/ICU, and day limit not reached)
                tmp[tmp[,prev_state]==ill_contag & tmp[,"d"] <sp$ill_dur & 
                    !tmp[,up_state]%in%c(ill_ICU,dead) ,up_state] <- ill_contag
                
                ## Update days (current count plus 1)
                tmp[tmp[,prev_state]==ill_contag & tmp[,up_state]==ill_contag,"d_new"] <- 
                    tmp[tmp[,prev_state]==ill_contag & tmp[,up_state]==ill_contag,"d"] + 1
                
                ## Getting healthy/immune if day limit reached
                tmp[tmp[,prev_state]==ill_contag & tmp[,"d"] >=sp$ill_dur & 
                    !tmp[,up_state]%in%c(ill_ICU,dead) ,up_state] <- immune
                
                ## Update days (reset to 1)
                tmp[tmp[,prev_state]==ill_contag & tmp[,up_state]==immune,"d_new"] <- 1

                ## #############################################################
                ## State: Ill, ICU #############################################
                
                ## Dying: Population at risk
                tmp$count   <- tmp[,prev_state]%in%c(ill_ICU)
                
                ##at_risk <- aggregate(count~dist_id,data=tmp,sum)
                at_risk <- group_by(tmp, dist_id) %>% summarise(count = sum(count))
                
                at_risk <- sum(at_risk$count)
                
                ## Check if anyone can die
                if(at_risk>0) {
                    
                    ## Potentially dying in ICU: Get their ages
                    age_sex <- tmp[tmp[,prev_state]%in%c(ill_ICU),c("age_gr","sex")]
                    age_sex <- paste(age_sex$age_gr,age_sex$sex)
                    
                    ## Match with ages for ICU risk of dying
                    age_sex <- match(age_sex,paste(surv_ICU$age,surv_ICU$sex))
                    
                    ## Assign individual risks
                    surv_ICU_i <- surv_ICU$surv_ICU[age_sex]
                    
                    ## Dying: Who dies? Sample
                    die <- as.numeric(runif(at_risk)>=surv_ICU_i)
                    
                    ## Recode
                    die[die == 0]  <- ill_contag
                    die[die == 1]  <- dead
                    ##die <- recode(die,'0' = ill_contag, '1' = dead)
                    
                    ## Dying: Assign state
                    tmp[tmp[,prev_state]==ill_ICU,up_state] <- die
                    
                }
                
                ## Staying ill/ICU (if not dead, and day limit not reached)
                tmp[tmp[,prev_state]==ill_ICU & tmp[,"d"] <lhc[it.ss,"icu_dur"] & 
                    !tmp[,up_state]%in%c(dead) ,up_state] <- ill_ICU
                
                ## Update days (current count plus 1)
                tmp[tmp[,prev_state]==ill_ICU & tmp[,up_state]==ill_ICU,"d_new"] <- 
                    tmp[tmp[,prev_state]==ill_ICU & tmp[,up_state]==ill_ICU,"d"] + 1
                
                ## Getting healthy/ICU if day limit reached
                tmp[tmp[,prev_state]==ill_ICU & tmp[,"d"] >=lhc[it.ss,"icu_dur"] & 
                    !tmp[,up_state]%in%c(dead) ,up_state] <- immune
                
                ## Update days (reset to 1)
                tmp[tmp[,prev_state]==ill_ICU & tmp[,up_state]==immune,"d_new"] <- 1
                
                ## If counter not increased (e.g., if dead)
                tmp[is.na(tmp[,"d_new"]),"d_new"] <- missing
                
                ## Move to main data frame
                sim[!sim[,prev_state]%in%c(dead,immune),c(up_state,"d")] <- tmp[c(up_state,"d_new")]
                
                ## Immune stay immune
                sim[sim[,prev_state]==immune,up_state] <- immune
                
                ## The dead stay dead
                sim[sim[,prev_state]==dead,up_state] <- dead

                ## #############################################################
                ## Calculate total daily numbers by counties ###################
                ## 0 ## Healthy
                cbtmp <- data.table(sim[sim[,up_state]==healthy,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                healthy_cases[,timestep] <- final_count
                
                ## 1 ## Infected, not contagious
                cbtmp <- data.table(sim[sim[,up_state]==inf_noncon,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                inf_noncon_cases[,timestep] <- final_count
                
                ## 2 ## Infected, contagious
                cbtmp <- data.table(sim[sim[,up_state]==inf_contag,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                inf_contag_cases[,timestep] <- final_count
                
                ## 3 ## Ill, contagious
                cbtmp <- data.table(sim[sim[,up_state]==ill_contag,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                ill_contag_cases[,timestep] <- final_count
                
                ## 4 ## Ill, ICU
                cbtmp <- data.table(sim[sim[,up_state]==ill_ICU,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                ill_ICU_cases[,timestep] <- final_count

                ## 5 ## Immune
                cbtmp <- data.table(sim[sim[,up_state]==immune,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                immune_cases[,timestep] <- final_count
                
                ## 6 ## Dead
                cbtmp <- data.table(sim[sim[,up_state]==dead,c(up_state,"dist_id")])
                case_count <- cbtmp[, .N, by=dist_id]
                final_count <- numeric(length(sim.struc$counties))
                final_count[match(case_count$dist_id,sim.struc$counties)] <- case_count$N
                dead_cases[,timestep] <- final_count

                ## #############################################################
                ## Flip-Flop sim timestep columns ##############################
                names(sim)[names(sim)==prev_state]  <- paste0("t",timestep+1)
                
                ## -----------------------------------------------------------------
                ## Checkpoint
                if ((timestep == ep$cp.time) & ep$cp.write) {
                    
                    cp.data <-  list(healthy_cases    = healthy_cases[,c(1:timestep)],
                                     inf_noncon_cases = inf_noncon_cases[,c(1:timestep)],
                                     inf_contag_cases = inf_contag_cases[,c(1:timestep)],
                                     ill_contag_cases = ill_contag_cases[,c(1:timestep)],
                                     ill_ICU_cases    = ill_ICU_cases[,c(1:timestep)],
                                     immune_cases     = immune_cases[,c(1:timestep)],
                                     dead_cases       = dead_cases[,c(1:timestep)],
                                     inf_cases        = inf_cases[,c(1:timestep)],
                                     icu_cases        = icu_cases[,c(1:timestep)],
                                     dead_cases_bICU  = dead_cases_bICU[,c(1:timestep)],
                                     mod_inf_cases    = mod_inf_cases[,c(1:timestep)],
                                     org_noncon_cases = org_noncon_cases[,c(1:timestep)],
                                     lockdowns        = lockdowns[,c(1:timestep)],
                                     sim              = sim)

                    saveRDS(cp.data,
                            paste0(ep$output,"/checkpoint.data-proc",it.ss,".RDS"))

                    print(paste("Checkpoint saved for up_state =",up_state,".",
                                "FlipFlop of sim columns is already done.",
                                "names(sim) =",paste(names(sim),collapse=",")))
                    
                }
                
            }
### <<< End of Simulation Model Loop ###########################################
################################################################################

            if (sp$results == "ALL") {
                ## -------------------------------------------------------------
                ## Results to put into the foreach result - fer                
                res.list <-  list(healthy_cases    = healthy_cases,   
                                  inf_noncon_cases = inf_noncon_cases,
                                  inf_contag_cases = inf_contag_cases,
                                  ill_contag_cases = ill_contag_cases,
                                  ill_ICU_cases    = ill_ICU_cases,
                                  immune_cases     = immune_cases,
                                  dead_cases       = dead_cases,
                                  inf_cases        = inf_cases,
                                  icu_cases        = icu_cases,
                                  dead_cases_bICU  = dead_cases_bICU,
                                  iter             = lhc[it.ss,],
                                  dist_id          = sim.struc$counties,
                                  mod_inf_cases    = mod_inf_cases,
                                  org_noncon_cases = org_noncon_cases)
            } else {
                res.list <-  list(ill_ICU_cases    = ill_ICU_cases,
                                  iter             = lhc[it.ss,],
                                  dist_id          = sim.struc$counties)
                                  
            }
            
            
            ## -----------------------------------------------------------------
            ## Return
            return(res.list)
            
        }
###### <<< End of Iteration over parameter space ###############################
################################################################################
    
    print("Parameter Iteration done")

################################################################################
### Export output and plot data if exec.procedure is "Basic-Param"            ##
################################################################################
    if (ep$exec.procedure == "Basic-Param") {
       
        ## Reformat foreach result fer to save a list of aggregated states per county
        if (sp$results == "ALL") {
            rr <- list(healthy=rbindlist(lapply(fer,
                                                function(x){data.frame(x$healthy_cases,
                                                                       x$iter,
                                                                       x$dist_id)})),
                       inf_noncon=rbindlist(lapply(fer,
                                                   function(x){data.frame(x$inf_noncon_cases,
                                                                          x$iter,
                                                                          x$dist_id)})),
                       inf_contag=rbindlist(lapply(fer,
                                                   function(x){data.frame(x$inf_contag_cases,
                                                                          x$iter,
                                                                          x$dist_id)})),
                       ill_contag=rbindlist(lapply(fer,
                                                   function(x){data.frame(x$ill_contag_cases,
                                                                          x$iter,
                                                                          x$dist_id)})),
                       ill_ICU=rbindlist(lapply(fer,
                                                function(x){data.frame(x$ill_ICU_cases,
                                                                       x$iter,
                                                                       x$dist_id)})),
                       immune_cases=rbindlist(lapply(fer,
                                                     function(x){data.frame(x$immune_cases,
                                                                            x$iter,
                                                                            x$dist_id)})),
                       dead_cases=rbindlist(lapply(fer,
                                                   function(x){data.frame(x$dead_cases,
                                                                          x$iter,
                                                                          x$dist_id)})),
                       inf_cases=rbindlist(lapply(fer,
                                                  function(x){data.frame(x$inf_cases,
                                                                         x$iter,
                                                                         x$dist_id)})),
                       icu_cases=rbindlist(lapply(fer,
                                                  function(x){data.frame(x$icu_cases,
                                                                         x$iter,
                                                                         x$dist_id)})),
                       dead_cases_bICU=rbindlist(lapply(fer,
                                                        function(x){data.frame(x$dead_cases_bICU,
                                                                               x$iter,
                                                                               x$dist_id)})),
                       mod_inf_cases=rbindlist(lapply(fer,
                                                      function(x){data.frame(x$mod_inf_cases,
                                                                             x$iter,
                                                                             x$dist_id)})),
                       org_noncon_cases=rbindlist(lapply(fer,
                                                         function(x){data.frame(x$org_noncon_cases,
                                                                                x$iter,
                                                                                x$dist_id)})))
        } else {
            
            rr <- list(ill_ICU=rbindlist(lapply(fer,
                                            function(x){data.frame(x$ill_ICU_cases,
                                                                   x$iter,
                                                                   x$dist_id)})))
        }
        
                   
        save(file=paste(ep$output.dir,
                        "Aggregated_Results-iter=",sp$iter,
                        "-sam_size=",pspace[["sam_size"]]$param[1],"-",
                        ep$export_name,".RData",sep=""),
             list="rr")

        
        ## #####################################################################
        ## Instant Plots                                                       #
        ## #####################################################################
        
        ## #####################################################################
        ## Cumulative plots of all states , Infections and ICU cases (total)
        ## across all parameter sets across the country
        if (sp$gplots) {
            outfile <- paste0(ep$output.dir,
                          "Rplots-ByCountry-totals-iter=",sp$iter,"-sam_size=",
                          pspace[["sam_size"]]$param[1],"-",
                          ep$export_name,".pdf",sep="")

        plots.by.country (outfile         = outfile,
                          sp              = sp,
                          seed_icu        = op$opt.target$icu,
                          seed_dea        = op$opt.target$dea,
                          iol             = iol,
                          pspace          = pspace,
                          rr              = rr,
                          ind.states   = c(1,2,3,4,5,6,7),
                          global.plot     = TRUE)
        }
        
        ## #####################################################################
        ## Cumulative plots of all states, aggregated once across the first
        ## column in the latin hypercube, across each direct parameter with
        ## more than one value and once across the parameter set of the
        ## directv parameter for R0effect.
        if (sp$cplots) {
            outfile <- paste0(ep$output.dir,
                              "Rplots-ByCountry-params-iter=",sp$iter,"-sam_size=",
                              pspace[["sam_size"]]$param[1],"-",
                              ep$export_name,".pdf",sep="")

            plots.by.country (outfile         = outfile,
                              sp              = sp,
                              seed_icu        = op$opt.target$icu,
                              seed_dea        = op$opt.target$dea,
                              iol             = iol,
                              pspace          = pspace,
                              rr              = rr,
                              ind.states   = c(1,2,3,4,5,6,7),
                              global.plot     = FALSE)    
        }
        
        ## #####################################################################
        ## Cumulative plots for every state, aggregated once across the
        ## first column in the latin hypercube, across each direct parameter
        ## with more than one value and once across the parameter set of the
        ## directv parameter for R0effect.
        if (sp$cplots.states) {
            
            ## Every diagram its own scale ------------------------------
            outfile <- paste0(ep$output.dir,
                              "Rplots-ByState-iter=",sp$iter,"-sam_size=",
                              pspace[["sam_size"]]$param[1],"-RMS-",
                              ep$export_name,".pdf",sep="")
            
            plots.by.state(outfile      = outfile,
                           sp           = sp,
                           seed_icu     = op$opt.target$icu.bs,
                           seed_dea     = op$opt.target$dea.bs,
                           iol          = iol,
                           pspace       = pspace,
                           rr           = rr,
                           region       = "state",
                           fix.lim      = FALSE,
                           filtered     = FALSE,
                           Sec.Axis=c("RMS"),
                           silent=FALSE,
                           fk.cases=rep(1/ 7,  7),
                           fk.sec  =rep(1/15, 15),
                           ind.states   = c(5))
        
            ## Every diagram its own scale with R0effects ----------------------
            ## This is only meaningfull if R0effects vary on state level -------
            if ( R0effect.region == "state" ) {
                outfile <- paste0(ep$output.dir,
                                  "Rplots-ByState-iter=",sp$iter,"-sam_size=",
                                  pspace[["sam_size"]]$param[1],"-R0effects-",
                                  ep$export_name,".pdf",sep="")
                
                plots.by.state(outfile      = outfile,
                               sp           = sp,
                               seed_icu     = iol$icu.cases.by.state,
                               seed_dea     = op$opt.target$dea.bs,
                               iol          = iol,
                               pspace       = pspace,
                               rr           = rr,
                               region       = "state",
                               fix.lim      = FALSE,
                               filtered     = FALSE,
                               Sec.Axis     = c("R0effect"),
                               silent       = FALSE,
                               fk.cases     = rep(1/ 7,  7),
                               fk.sec       = rep(1/15, 15),
                               ind.states   = c(5))
            }
            
            ## Global fixed scale accross all diagramms -----------------
            ## outfile <- paste0(ep$output.dir,
            ##                   "Rplots-ByState-iter=",sp$iter,"-sam_size=",
            ##                   pspace[["sam_size"]]$param[1],"-",
            ##                   ep$export_name,"-fixed-scale",".pdf",sep="")
            ## 
            ## plots.by.state(outfile      = outfile,
            ##                sp           = sp,
            ##                seed_icu     = op$opt.target$icu.bs,
            ##                seed_dea     = op$opt.target$dea.bs,
            ##                iol          = iol,
            ##                pspace       = pspace,
            ##                rr           = rr,
            ##                region       = "state",
            ##                fix.lim      = TRUE,
            ##                ind.states   = c(1,2,3,4,5,6,7))
        }
        
        ## #####################################################################
        ## Cumulative plots for every state, aggregated once across the
        ## first column in the latin hypercube, across each direct parameter
        ## with more than one value and once across the parameter set of the
        ## directv parameter for R0effect.
        if (sp$cplots.nuts2) {

            ## Every diagram its own scale ------------------------------
            outfile <- paste0(ep$output.dir,
                              "Rplots-ByNUTS2-iter=",sp$iter,"-sam_size=",
                              pspace[["sam_size"]]$param[1],"-RMS-",
                              ep$export_name,".pdf",sep="")
            
            plots.by.state(outfile      = outfile,
                           sp           = sp,
                           seed_icu     = op$opt.target$icu.nuts2,
                           seed_dea     = NULL,
                           iol          = iol,
                           pspace       = pspace,
                           rr           = rr,
                           region       = "nuts2",
                           fix.lim      = FALSE,
                           filtered = FALSE,
                           Sec.Axis=c("RMS"),
                           silent=FALSE,
                           ind.states   = c(1,2,3,4,5,6,7))
            
            ## Every diagram its own scale with R0effects ---------------
            ## This is only meaningfull if R0effects vary on nuts2 level -------
            if ( R0effect.region == "nuts2" ) {
                
                outfile <- paste0(ep$output.dir,
                                  "Rplots-ByNUTS2-iter=",sp$iter,"-sam_size=",
                                  pspace[["sam_size"]]$param[1],"-R0effect-",
                                  ep$export_name,".pdf",sep="")
                
                plots.by.state(outfile      = outfile,
                               sp           = sp,
                               seed_icu     = op$opt.target$icu.nuts2,
                               seed_dea     = NULL,
                               iol          = iol,
                               pspace       = pspace,
                               rr           = rr,
                               region       = "nuts2",
                               fix.lim      = FALSE,
                               filtered = FALSE,
                               Sec.Axis=c("R0effect"),
                               ind.states   = c(1,2,3,4,5,6,7))
            }
            
            ## Global fixed scale accross all diagramms -----------------
            outfile <- paste0(ep$output.dir,
                              "Rplots-ByNUTS2-iter=",sp$iter,"-sam_size=",
                              pspace[["sam_size"]]$param[1],"-fixed-scale-",
                              ep$export_name,".pdf",sep="")
            
            plots.by.state(outfile      = outfile,
                           sp           = sp,
                           seed_icu     = op$opt.target$icu.nuts2,
                           seed_dea     = NULL,
                           iol          = iol,
                           pspace       = pspace,
                           rr           = rr,
                           region       = "nuts2",
                           fix.lim      = TRUE,
                           ind.states   = c(1,2,3,4,5,6,7))
        }
                
################################################################################
### If exec.procedure == "Optimization" calculate fittness function value   ####
################################################################################
    } else { 
        
        ## Simulation timeframe ------------------------------------------------
        time <- seq(from=seed_date,to=seed_date+sp$time_n-1,by="days")

        if ( op$opt.target.region == "country" ) {
            ## -----------------------------------------------------------------
            ## Global target function ------------------------------------------
            
            if ( op$opt.target.deaths ) {
                
                ## Aggregate deaths by dist_id ---------------------------------
                col.Sums  <- lapply(fer,function(x){colSums(x[["dead_cases"]])})
                data.it <- as.data.frame(do.call(rbind, col.Sums))
                
                ### Calculate mean across the aggregated data ------------------
                df.gg <- data.frame(time=time, MeanOfAggregated = apply(data.it, 2, mean))
                
                ##--------------------------------------------------------------
                ## Daywise distance against the optimization target ------------
                diff.d <-
                    op$opt.target$dea$deaths[op$opt.target$dea$date %in% df.gg$time] -
                    df.gg                [df.gg$time %in% op$opt.target$dea$date,"MeanOfAggregated"]
                
                sddq  <- sqrt(sum(diff.d*diff.d)/length(diff.d))
                
            } else {
                sddq  <- 0
            }
            
            if ( op$opt.target.icu ) {
                
                ## Aggregate icu_cases by dist_id ------------------------------
                col.Sums  <- lapply(fer,function(x){colSums(x[["ill_ICU_cases"]])})
                data.it <- as.data.frame(do.call(rbind, col.Sums))
                
                ## Calculate mean across the aggregated data -------------------
                df.gg <- data.frame(time=time, MeanOfAggregated = apply(data.it, 2, mean))
                
                ##--------------------------------------------------------------
                ## Daywise distance against the optimization target ------------
                diff.i <- op$opt.target$icu$cases[op$opt.target$icu$date %in% df.gg$time] -
                    df.gg[df.gg$time %in% op$opt.target$icu$date,"MeanOfAggregated"]
                
                sdiq  <- sqrt(sum(diff.i*diff.i)/length(diff.i))
                
            } else {
                sdiq  <- 0
            }

        } else if ( op$opt.target.region == "state" ) {
            ## -----------------------------------------------------------------
            ## State local target function -------------------------------------

            if ( op$opt.target.deaths ) {

                sddq.l <- list()
                st.id <- unique(as.integer(fer[[1]]$dist_id/1000))
                
                for ( ii in seq(length(st.id)) ) {

                    ## Aggregate deaths by dist_id -----------------------------
                    

                    sel.c <- as.integer(fer[[1]]$dist_id/1000)==st.id[ii]
                    if ( sum( sel.c ) > 1) {
                        ## For States with more than one county ----------------
                        col.Sums  <- lapply(fer,function(x){
                            colSums(x$dead_cases[as.integer(x$dist_id/1000)==st.id[ii],])})
                    } else {
                        ## For States with only one county ---------------------
                        col.Sums  <- lapply(fer,function(x){
                            x$dead_cases[as.integer(x$dist_id/1000)==st.id[ii],]})
                    }
                    
                    data.it <- as.data.frame(do.call(rbind, col.Sums))
                    
                    ## Calculate mean across the aggregated data ---------------
                    df.gg <- data.frame(time=time, MeanOfAggregated = apply(data.it, 2, mean))
                    
                    ## Select target deaths by state ---------------------------
                    ot <- op$opt.target$dea.bs[op$opt.target$dea.bs$Bundesland==iol[["states"]][st.id[ii],"Name"],]
                    ot$date <- as.Date(ot$date)
                    
                    diff.d  <- ot   [ot$date %in% df.gg$time,"SummeTodesfaelle"] -
                               df.gg[df.gg$time %in% ot$date,"MeanOfAggregated"]
                    
                    sddq.l[[ii]] <- diff.d
                }

                sddq.ul <- unlist(sddq.l)
                sddq <- sqrt(sum(sddq.ul*sddq.ul)/length(sddq.ul))
                
            } else {
                sddq  <- 0
            }

            ## target Local ICU-Cases ------------------------------------------
            if ( op$opt.target.icu ) {
                
                sdiq.l <- list()
                st.id <- unique(as.integer(fer[[1]]$dist_id/1000))
                
                for ( ii in seq(length(st.id)) ) {
                    
                    ## Aggregate ICU cases by state via dist_id ----------------

                    ## Select counties of current state ------------------------
                    sel.c <- as.integer(fer[[1]]$dist_id/1000)==st.id[ii]

                    if ( sum( sel.c ) > 1) {
                        ## In case the state has  more than one county ---------
                        col.Sums  <- lapply(fer,function(x){
                            colSums(x$ill_ICU_cases[sel.c,])})

                    } else {
                        ## In case the state has only one county ---------------
                        col.Sums  <- lapply(fer,function(x){
                            x$ill_ICU_cases[sel.c,]})
                    }
                    data.it <- as.data.frame(do.call(rbind, col.Sums))
                
                    ## Calculate mean across the aggregated data ---
                    df.gg <- data.frame(time=time, MeanOfAggregated = apply(data.it, 2, mean))
                    
                    ## Select target icu cases by state -----------
                    ot <- op$opt.target$icu.bs[op$opt.target$icu.bs$state ==
                                               iol[["states"]][st.id[ii],"Shortcut"],]
                   
                    if (op$opt.filter) {
                        ## Filter Kernel --------
                        f7 <- rep(1/7, 7)
                        
                        ## Save original to replace NAs -----
                        f.cases  <- ot$cases
                        ## Filter ------------------------------------
                        ot$cases <- stats::filter(ot$cases,f7,sides=2)
                        ## Substitute NAs with original data ---------
                        w.na <- which(is.na(ot$cases))
                        ot$cases[w.na] <- f.cases[w.na]

                        ## Save original to replace NAs -----
                        f.cases  <- df.gg$MeanOfAggregated
                        ## Filter ------------------------------------
                        df.gg$MeanOfAggregated <- stats::filter(df.gg$MeanOfAggregated,
                                                                f7,sides=2)
                        ## Substitute NAs with original data ---------
                        w.na <- which(is.na(df.gg$MeanOfAggregated))
                        df.gg$MeanOfAggregated[w.na] <- f.cases[w.na]
                    }

                    diff.i  <- ot   [ot$date %in% df.gg$time,"cases"] -
                               df.gg[df.gg$time %in% ot$date,"MeanOfAggregated"]
                    
                    sdiq.l[[ii]] <- diff.i
                }
                
                sdiq.ul <- unlist(sdiq.l)
                sdiq <- sqrt(sum(sdiq.ul*sdiq.ul)/length(sdiq.ul))
                
            } else {
                sdiq  <- 0
            }
            
        } else if (op$opt.target.region == "nuts2") {

            ## target Local ICU-Cases ------------------------------------------
            if ( op$opt.target.icu ) {

                sdiq.l <- list()
                nuts2.id <- unique(iol$counties[iol$counties$dist_id %in%
                                                fer[[1]]$dist_id,"Nuts2"])
                
                for ( ii in seq(length(nuts2.id)) ) {
                    
                    ## Aggregate ICU cases by NUTS2 region via dist_id ---------
                    
                    ## Select counties of current NUTS2 region -----------------
                    sel.c <- iol$counties[iol$counties$dist_id %in%
                                          fer[[1]]$dist_id,"Nuts2" ] == nuts2.id[ii]

                    if ( sum( sel.c ) > 1) {
                        ## In case the NUTS2 region has  more than one county --
                        col.Sums  <- lapply(fer,function(x){
                            colSums(x$ill_ICU_cases[sel.c,])})
                        
                    } else {
                        ## In case the NUTS2 region has only one county --------
                        col.Sums  <- lapply(fer,function(x){
                            x$ill_ICU_cases[sel.c,]})
                    }
                    data.it <- as.data.frame(do.call(rbind, col.Sums))
                
                    ## Calculate mean across the aggregated data ---
                    df.gg <- data.frame(time=time, MeanOfAggregated = apply(data.it, 2, mean))

                    ## Select target icu cases by NUTS2 region ----
                    ot <- op$opt.target$icu.nuts2[op$opt.target$icu.nuts2$Nuts2 ==
                                                  nuts2.id[ii],]

                    if (op$opt.filter) {
                        ## Filter Kernel --------
                        f7 <- rep(1/7, 7)
                        
                        ## Save original to replace NAs -----
                        f.cases  <- ot$cases
                        ## Filter ------------------------------------
                        ot$cases <- stats::filter(ot$cases,f7,sides=2)
                        ## Substitute NAs with original data ---------
                        w.na <- which(is.na(ot$cases))
                        ot$cases[w.na] <- f.cases[w.na]

                        ## Save original to replace NAs -----
                        f.cases  <- df.gg$MeanOfAggregated
                        ## Filter ------------------------------------
                        df.gg$MeanOfAggregated <- stats::filter(df.gg$MeanOfAggregated,
                                                                f7,sides=2)
                        ## Substitute NAs with original data ---------
                        w.na <- which(is.na(df.gg$MeanOfAggregated))
                        df.gg$MeanOfAggregated[w.na] <- f.cases[w.na]
                                        }
                    
                    diff.i  <- ot   [ot$date %in% df.gg$time,"cases"] -
                               df.gg[df.gg$time %in% ot$date,"MeanOfAggregated"]
                    
                    sdiq.l[[ii]] <- diff.i
                }
               
                sdiq.ul <- unlist(sdiq.l)
                sdiq <- sqrt(sum(sdiq.ul*sdiq.ul)/length(sdiq.ul))

            } else {
                sdiq  <- 0
            }
        }

        if ((ep$exec.procedure == "Optimization") &
            (ep$parallel.method == "PSOCK")         ) {
            print("Closing Cluster")
            stopCluster(cls)
        ##    ##rm(list="cls")
        ##    ##gc()
        }                

        ## Return the sum of the quaratic daily distances ---
        return(-(sddq+sdiq))

    }
    ## -------------------------------------------------------------------------
    ## <<< End of Switch for Results Export ------------------------------------
    ## -------------------------------------------------------------------------

}
################################################################################
### <<< End of Simulation_Model Function                                    ####
################################################################################

################################################################################
#' CoSMic model function wrapper
#'
#' This function wraps the CoSMic model function so that it can be used in the
#' GA algorithm as the objective function.
#'
#' @param x Numeric vector with model parameters subject to optimization.
#' @param ep Execution parameter list. Use [set.exec.params()] in order to
#'           create a valid layout.
#' @param sp List with static model parameters. Use [set.static.params()] to 
#'           create a valid layout.
#' @param iol Input data list. Use [load.input()] to load needed fies and
#'            [init.connectivity()] in order to create a valid date layout.
#' @param pspace List holding the parameter space with potentially variable
#'               model parameters. Use the setter function [set.pspace()] to 
#'               add parameters.
#' @param sim.struc List with population data. Use [init.spatial.population()]
#'                  in order to create a valid layout.
#' @param op List with steering parameters for the optimization process. \cr
#'           Use [set.optimization.params()] in order to create a valid layout 
#'           and\cr [init.reference.data()] in order to init the optimization 
#'           targets based on observed data.
#'
#' @return A scalar value calculated according to the settings given in op.
#' 
ff <- function(x,
               ep, sp, iol, pspace, sim.struc, op) {

    a <- CoSMic (ep = ep,
                 sp = sp,
                 iol = iol,
                 pspace=pspace,
                 sim.struc = sim.struc,
                 op  = op,
                 opt = x)
        print(a)
    return(a)
}


################################################################################
#' GA algorithm monitoring function
#'
#' The function provides intermediate output after each iteration of the GA
#' algorithm.
#'
#' @param obj An object provided by the GA function.
#' @param digits The number of digits provided by \code{getOption("digits")}.
#' @param sp.int List with static model parameters as created by
#'               [set.static.params()].
#' @param op.int List with steering parameters for the optimization process
#'               as created by \cr [set.static.params()].
#' 
GA.Monitor <- function(obj, digits = getOption("digits"),
                       sp.int=static.params, op.int=opt.params) 
{
    print(paste("GA@iter:",obj@iter,sep=""))

    mon.file <- paste0("Best-Sol-iter=",sp.int$iter,"-sam_size=",
                      pspace[["sam_size"]]$param[1],"-",
                      sp.int$export_name,".csv",sep="")
    
    if( ! file.exists(mon.file) ) {
        
        write.table(t(op.int$opt.names),
                    mon.file,
                    sep=",", row.names=F, col.names=F)
    }
    
    if ( obj@iter > 1 ) {

        write.table(x = paste(obj@iter,
                              paste(obj@bestSol[[obj@iter-1]],sep="",collapse=",")),
                    file = mon.file,
                    sep=",", row.names=F, col.names=F, dec=".", append=T)
        
    }
}


################################################################################
#' Application of the GA algorithm to CoSMic.
#'
#' The function applies the GA algorithm to the CoSMic simulation model
#' function. it uses the wrapper function [ff()] as the objective function and
#' [GA.Monitor()] to return intermediate results during the course of the
#' optimization.
#'
#' @param ep Execution parameter list. Use [set.exec.params()] in order to
#'           create a valid layout.
#' @param sp List with static model parameters. Use [set.static.params()] to 
#'           create a valid layout.
#' @param iol Input data list. Use [load.input()] to load needed fies and
#'            [init.connectivity()] in order to create a valid date layout.
#' @param pspace List holding the parameter space with potentially variable
#'               model parameters. Use the setter function [set.pspace()] to 
#'               add parameters.
#' @param sim.struc List with population data. Use [init.spatial.population()]
#'                  in order to create a valid layout.
#' @param op List with steering parameters for the optimization process. \cr
#'           Use [set.optimization.params()] in order to create a valid layout 
#'           and\cr [init.reference.data()] in order to init the optimization 
#'           targets based on observed data.
#' @param cl A parallel cluster prepared by [init.parallel.execution()].
#'
#' @export
CoSMic.Opt <- function(ep, sp, iol, pspace, sim.struc, op, cl) {
    
    ## Execute Genetic Algorithm --------------------------------------
    registerDoRNG(as.integer(Sys.time()))
    
    GA<-ga(type="real-valued",
           fitness=ff,
           ep = ep,
           sp = sp,
           iol = iol,
           pspace=pspace,
           sim.struc = sim.struc,
           op=op,
           lower=op$opt.lb,
           upper=op$opt.ub,
           popSize= op$opt.pop.size,
           maxiter= op$opt.max.iter,
           monitor = GA.Monitor,
           keepBest= TRUE,
           suggestions = op$use.sug.sol,
           parallel=cl)        
    
    ##--------------------------------------------------------------------------
    ## Evaluation and Results Export -------------------------------------------
    ##--------------------------------------------------------------------------

    ## Show, and save optimization result ----------------------------
    print(summary(GA))
      
    save(GA,file=paste("GA-iter=",sp$iter,"-sam_size=",pspace[["sam_size"]]$param[1],
                       "-",ep$export_name,".RData",sep=""))
    
}

################################################################################
#' Prepare parameter space
#'
#' The function initializes the data.frame carrying the different sets of model
#' parameters resulting from the parameter variations set in the pspace list.
#' 
#' @param pspace The parameter list pspace set by repeated calls to [set.pspace()]
#' @param A list with static model parameters as described in [set.static.params()].
#'
#' @return A data.frame with dimension
#'         \[<# different evaluations> x <potentially_variable_model_params>\]
#'         If all model parameters in pspace are fixed, i.e. not variable
#'         dim(lhc) will be \[sp$iter x <potentially_variable_model_params>\]
#' 
#' @export
init.lhc <-  function(pspace,sp) {
        
    ## Extract parameter types from paramter list ---------------------
    pspace.types <- lapply(pspace,function(x){x$type})
    
    # Extract numeric type of parameters -----------------------------
    pspace.num.types <- lapply(pspace,function(x){x$num.type})
        
    # Number of parameters of type direct ----------------------------
    n.direct  <- sum(pspace.types=="direct")
    
    # Number of parameters of type direct ----------------------------
    n.directv  <- sum(pspace.types=="directv")

    # Number of parameters of type direct ----------------------------
    n.directl  <- sum(pspace.types=="directl")
    
    # Extract parameters of type dist --------------------------------
    pspace.params <- lapply(pspace[which(pspace.types %in% c("dist","distv"))],
                            function(x){x$param})
        
    # Extract standard deviations of type dist -----------------------
    pspace.sd <- lapply(pspace[which(pspace.types %in% c("dist","distv"))],
                        function(x){x$sd})
    
    # Number of elements to vary in parameters of type dist ----------
    n.dist    <-  sum(unlist(pspace.sd)>0)
    
    # Calculate Random latin hyper cube for all paramters with -------
    # type = "dist"
    if (n.dist > 0) {
        
        lhc <- as.data.frame(randomLHS(n=sp$lhc.samples, k=n.dist))
        
        ## Scale columns of lhc to param mean +- sd ------------------
        lhc <- as.data.frame(mapply(
            function(x,y,z){ x + (y*2-1) * z/100 * x },
            ## x = All parameters of type dist with sd > 0 ---
            ##     Is reused by mapply for all rows in lhc ---
            unlist(pspace.params)[which(unlist(pspace.sd)>0)],
            ## each row of lhc ---
            lhc,
            ## z = All sd from parameters of type dist with sd > 0 ---
            ##     Is reused by mapply for all rows in lhc         ---
            unlist(pspace.sd)[which(unlist(pspace.sd)>0)]))
            
        ## Add colums to lhc for parameters of type dist with sd = 0 -----------
        ## unlist(pspace.params)[which(unlist(pspace.sd)==0)]) returns all
        ## parameters of type dist with sd = 0
        ## rep(as.data.frame( ... ,dim(lhc)[1])) returns a list of one column 
        ## data.frames
        ## do.call(rbind, ... ) binds the elements of the rep list as rows
        ## of a 2D numeric.
        lhc <- data.frame(lhc,    
                          do.call(rbind,
                                  rep(as.data.frame(
                                      unlist(pspace.params)[which(unlist(pspace.sd)==0)]
                                  ),
                                  dim(lhc)[1])
                                  ))
        names(lhc) <- c(
            names( unlist(pspace.params[[1]])[which(unlist(pspace.sd) >0)]),
            names( unlist(pspace.params[[1]])[which(unlist(pspace.sd)==0)]) )
        
        rownames(lhc)  <- NULL

        ## Add 0 to week ids less than 10 --------------------
        if (sum(is.na(as.numeric(substring(names(lhc),5)))) == 0) {
            names(lhc) <- apply(
                data.frame(
                    substring(names(lhc),1,4),
                    sprintf("%02d", as.numeric(substring(names(lhc),5)))
                ),1,paste,collapse="")
        }
        
        ## Order the columns --------------------------------
        if (dim(lhc)[2] > 1) {
            lhc <- lhc[,order(names(lhc))]
        }
        
        ## randomLHS returns floating point values. ----------------------------
        ## Convert values with num.type = "int"
        ## int.params <- names(pspace.num.types)[which(pspace.num.types=="int")]
        ## lhc[,int.params] <- round(lhc[,int.params],digits=0)
        
        if ( n.direct > 0) {
            dp <- pspace[[which(pspace.types=="direct")[1]]]$param
            lhc <- data.frame(
                lhc %>% slice(rep(row_number(), length(dp))),
                rep(dp,each=dim(lhc)[1]))
        }
        
    } else {
            
        ## Init lhc ------------------------------------------------------------
            
        if ( sum(which(pspace.types=="dist")) > 0 ) {
            ## If "dist" type parameters are given but all have sd = 0 ---------

            lhc <- data.frame( t(unlist(lapply(pspace[which(pspace.types=="dist")],
                                               function(x){x$param}))),
                              pspace[[which(pspace.types=="direct")[1]]]$param)
            
        } else {
            ## If no "dist" type parameters are given --------------------------
            lhc <- data.frame(pspace[[which(pspace.types=="direct")[1]]]$param)
        }
    }
    
    ## Direct scalar parameters to lhc -----------------------------------------
    if ( n.direct > 0 ) {
        
        ## Permutate the random latin hyper cube for the paramters ---
        ## with type = "direct"
        for ( i in pspace[which(pspace.types=="direct")[-1]] ) {
            
            lhc <- data.frame(lhc %>% slice(rep(row_number(), length(i$param))),
                              rep(i$param,each=dim(lhc)[1]))
        }
    }
    
    ## Names of direct parameters ------------------------------------
    names(lhc)[seq(dim(lhc)[2]-n.direct+1,dim(lhc)[2])] <-
        names(pspace.types[which(pspace.types=="direct")])
    
    ## Direct vector parameters to lhc -----------------------------------------
    if ( n.directv > 0 ) {
        
        ## Permutate the random latin hyper cube for the paramters  --
        ## with type = "directv"
        for ( i in pspace[which(pspace.types=="directv")] ) {
            lhc <- data.frame(lhc %>% slice(rep(row_number(), each=dim(i$param)[1])),
                              i$param)
        }
    }
    
    ## Direct list parameters to lhc -------------------------------------------
    if ( n.directl > 0 ) {
        
        ## Permutate the random latin hyper cube for the paramters ---
        ## with type = "directl"
        for ( i in pspace[which(pspace.types=="directl")] ) {
            ## Convert the list of data.frames to a data.frame with
            ## as many rows as data.frames are in the list and as
            ## many columns as each data.frame has elements
            lhc.dat <- do.call(rbind,lapply(i$param,unlist))
            
            rownames(lhc.dat) <- NULL
            
            lhc <- data.frame(lhc %>% slice(rep(row_number(), each=dim(lhc.dat)[1])),
                              lhc.dat)
        }
    }
    
    ## Each Parameter set sp$iter times ---------------------------------
    lhc <- lhc %>% slice(rep(row_number(), sp$iter))
    
    ## Order according to equal parameter sets -----------------------
    lhc <- arrange_all(lhc)
    
    ## Add sp$iter as a parameter set column ----------------------------
    lhc <- data.frame(lhc,iter=rep(1:sp$iter,length.out=dim(lhc)[1]))

    return(lhc)
}
