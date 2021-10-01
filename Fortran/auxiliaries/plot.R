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
src.dir <- "/quobyte/qb1/s33094/s33094/cosmic/CoSMic/R"
for ( i in dir(src.dir,pattern="*.R$") ) { source(paste(src.dir,i,sep="/")) }

args = commandArgs(trailingOnly=TRUE)

date.str  <- args[1]
input.dir <- args[2]

cur.dir <- getwd()

################################################################################
### Reload baseline parameters                                                 #
################################################################################
sp            <- readRDS(paste(input.dir,dir(input.dir,pattern="static"),sep="/"))
iol           <- readRDS(paste(input.dir,dir(input.dir,pattern="input"),sep="/"))
pspace        <- readRDS(paste(input.dir,dir(input.dir,pattern="pspace"),sep="/"))

################################################################################
### Load Fortran results                                                       #
################################################################################
setwd("../output-1P/")

icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")
dead       <- read.table(paste0(date.str,"dead_cases.csv"),head=TRUE,sep=",")
immune     <- read.table(paste0(date.str,"immune_cases.csv"),head=TRUE,sep=",")

rr<-list("healthy"=healthy,"inf_noncon"=inf_noncon,
         "inf_contag"=inf_contag,"ill_contag"=ill_contag,
         "ill_ICU"=icu,"dead"=dead,"immune"=immune)

setwd("../output-10p/")
icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")
dead       <- read.table(paste0(date.str,"dead_cases.csv"),head=TRUE,sep=",")
immune     <- read.table(paste0(date.str,"immune_cases.csv"),head=TRUE,sep=",")

rr <- list("healthy"=rbind(rr$healthy,healthy),
           "inf_noncon"=rbind(rr$inf_noncon,inf_noncon),
           "inf_contag"=rbind(rr$inf_contag,inf_contag),
           "ill_contag"=rbind(rr$ill_contag,ill_contag),
           "ill_ICU"=rbind(rr$ill_ICU,icu),
           "dead"=rbind(rr$dead,dead),
           "immune"=rbind(rr$immune,immune))


setwd("../output-20p/")
icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")
dead       <- read.table(paste0(date.str,"dead_cases.csv"),head=TRUE,sep=",")
immune     <- read.table(paste0(date.str,"immune_cases.csv"),head=TRUE,sep=",")

rr <- list("healthy"=rbind(rr$healthy,healthy),
           "inf_noncon"=rbind(rr$inf_noncon,inf_noncon),
           "inf_contag"=rbind(rr$inf_contag,inf_contag),
           "ill_contag"=rbind(rr$ill_contag,ill_contag),
           "ill_ICU"=rbind(rr$ill_ICU,icu),
           "dead"=rbind(rr$dead,dead),
           "immune"=rbind(rr$immune,immune))


setwd("../output-50p/")
icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")
dead       <- read.table(paste0(date.str,"dead_cases.csv"),head=TRUE,sep=",")
immune     <- read.table(paste0(date.str,"immune_cases.csv"),head=TRUE,sep=",")

rr <- list("healthy"=rbind(rr$healthy,healthy),
           "inf_noncon"=rbind(rr$inf_noncon,inf_noncon),
           "inf_contag"=rbind(rr$inf_contag,inf_contag),
           "ill_contag"=rbind(rr$ill_contag,ill_contag),
           "ill_ICU"=rbind(rr$ill_ICU,icu),
           "dead"=rbind(rr$dead,dead),
           "immune"=rbind(rr$immune,immune))

setwd("../output-100p/")
icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")
dead       <- read.table(paste0(date.str,"dead_cases.csv"),head=TRUE,sep=",")
immune     <- read.table(paste0(date.str,"immune_cases.csv"),head=TRUE,sep=",")

rr <- list("healthy"=rbind(rr$healthy,healthy),
           "inf_noncon"=rbind(rr$inf_noncon,inf_noncon),
           "inf_contag"=rbind(rr$inf_contag,inf_contag),
           "ill_contag"=rbind(rr$ill_contag,ill_contag),
           "ill_ICU"=rbind(rr$ill_ICU,icu),
           "dead"=rbind(rr$dead,dead),
           "immune"=rbind(rr$immune,immune))

plots.by.country (outfile         = "plot.by.country.pdf",
                  sp              = sp,
                  seed_icu        = iol$icu.cases.by.country,
                  seed_dea        = iol$dead.cases.by.country,
                  iol             = iol,
                  pspace          = pspace,
                  rr              = rr,
                  ind.states      = c(1,2,3,4,5,7,6),
                  global.plot     = FALSE,
                  split.in = "sam_size")

## Every diagram its own scale ------------------------------
plots.by.state(outfile      = "plot.by.state.pdf",
               sp           = sp,
               seed_icu     = iol$icu.cases.by.state,
               seed_dea     = NULL,
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
               ind.states   = c(1,2,3,4,5,7,6),
               split.in = "sam_size")
