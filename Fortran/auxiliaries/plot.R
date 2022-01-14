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
################################################################################
################################################################################

## Better safe than sorry ------------------------------------------------------
rm(list=ls())

### Libraries ##################################################################

library(dplyr)
## Switch off dplyr grouping message ----------------------
## See https://dplyr.tidyverse.org/reference/summarise.html
options(dplyr.summarise.inform=FALSE) # 

library(rlist)
library(data.table)

library(ggplot2)    # Plotting                 ---
library(lhs)        # Latin Hypercube Sampling ---
library(GA)         # Genetic Algorithms       ---
library(doRNG)      # Parallel random seed     ---
library(grid)       # Gridplots                ---
library(gridExtra)
library(pracma)
library(RColorBrewer)

### Import functions ###########################################################
base.dir  <- ""
src.dir   <- paste0(base.dir,"/quobyte/qb1/s33094/s33094/cosmic/CoSMic/R")
for ( i in dir(src.dir,pattern="*.R$") ) { source(paste(src.dir,i,sep="/")) }

#args = commandArgs(trailingOnly=TRUE)
#date.str  <- args[1]
#input.dir <- args[2]
input.dir <- dir(pattern="Results")

################################################################################
### Reload baseline parameters                                                 #
################################################################################
sp            <- readRDS(paste(input.dir,dir(input.dir,pattern="static"),sep="/"))
iol           <- readRDS(paste(input.dir,dir(input.dir,pattern="input"),sep="/"))
pspace        <- readRDS(paste(input.dir,dir(input.dir,pattern="pspace"),sep="/"))

healthy    <- fres.to.dataframe(input.dir,"healthy_cases_")
inf_noncon <- fres.to.dataframe(input.dir,"inf_noncon_cases_")
inf_contag <- fres.to.dataframe(input.dir,"inf_contag_cases_")
ill_contag <- fres.to.dataframe(input.dir,"ill_contag_cases_")
ill_ICU    <- fres.to.dataframe(input.dir,"ill_ICU_cases_")
dead       <- fres.to.dataframe(input.dir,"dead_cases_")
immune     <- fres.to.dataframe(input.dir,"immune_cases_")
                                                                
rr<-list(healthy=healthy,       inf_noncon=inf_noncon,
         inf_contag=inf_contag, ill_contag=ill_contag,
         ill_ICU=ill_ICU,       dead=dead,
         immune=immune)

plots.by.country (outfile         = "plot.by.country.pdf",
                  sp              = sp,
                  seed_icu        = iol$icu.cases.by.country,
                  seed_dea        = iol$dead.cases.by.country,
                  iol             = iol,
                  pspace          = pspace,
                  rr              = rr,
                  ind.states      = NULL,
                  global.plot     = FALSE,
                  split.in        ="SH1")

## Every diagram its own scale ------------------------------
plots.by.state(outfile      = "plot.by.state.pdf",
               sp           = sp,
               seed_icu     = iol$icu.cases.by.state,
               seed_dea     = iol$dead.cases.by.state,
               iol          = iol,
               pspace       = pspace,
               rr           = rr,
               region       = "state",
               fix.lim      = FALSE,
               filtered     = FALSE,
               Sec.Axis     = NULL,
               silent       = FALSE,
               ind.states   = c(5),
               relative     = FALSE,
               split.in     = "SH1",
               prog         = NULL)
