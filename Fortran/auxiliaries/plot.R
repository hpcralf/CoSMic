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
base.dir  <- "/run/user/1000/gvfs/sftp:host=vulcan.hww.hlrs.de"
src.dir   <- paste0(base.dir,"/quobyte/qb1/s33094/s33094/cosmic/CoSMic/R")
input.dir <- paste0(base.dir,"/lustre/nec/ws3/ws/hpcralf-BIB/cosmic-projection-W64+6-06-02-21-20:07:15/Results-v12.0-2021-06-02_20:07:17")
data.dir  <- "./"

for ( i in dir(src.dir,pattern="*.R$") ) { source(paste(src.dir,i,sep="/")) }

#args = commandArgs(trailingOnly=TRUE)
#date.str  <- args[1]
#input.dir <- args[2]

################################################################################
### Reload baseline parameters                                                 #
################################################################################
sp            <- readRDS(paste(input.dir,dir(input.dir,pattern="static"),sep="/"))
iol           <- readRDS(paste(input.dir,dir(input.dir,pattern="input"),sep="/"))
pspace        <- readRDS(paste(input.dir,dir(input.dir,pattern="pspace"),sep="/"))

################################################################################
### Load Fortran results                                                       #
################################################################################
ldf <- list()
cnt <- 0
for ( f in  dir("./",pattern="ill_ICU_cases.csv")) {

    to.read = file(f,"rb")

    vec<-readBin(to.read, integer(), size=4, n=401*106*24*2, endian = "little")

    mat<-matrix(vec[1:(401*106*24*2)],nrow=401)

    df <- data.frame(mat)

    ldf[[f]] <- do.call(rbind,apply(do.call(cbind,df),1,function(x){data.frame(t(matrix(x,nrow=106)))}))
    
    ldf[[f]]$SH1       <- rep(rep(c(1,2),each=24),401)+cnt
    ldf[[f]]$iter      <- rep(seq(24),401*2)
    ldf[[f]]$x.dist_id <- rep(iol$counties$dist_id,each=24*2)

    ldf[[f]]<-ldf[[f]][order(ldf[[f]]$SH1,ldf[[f]]$iter),]
    cnt <- cnt + 2
}

tmp<-do.call(rbind,ldf)
tmp<-tmp[order(tmp$SH1,tmp$iter),]

rr<-list("healthy"=tmp,"inf_noncon"=2,
         "inf_contag"=3,"ill_contag"=4,
         "ill_ICU"=tmp,"dead"=6,"immune"=7)

plots.by.country (outfile         = "plot.by.country.pdf",
                  sp              = sp,
                  seed_icu        = iol$icu.cases.by.country,
                  seed_dea        = iol$dead.cases.by.country,
                  iol             = iol,
                  pspace          = pspace,
                  rr              = rr,
                  ind.states      = c(5),
                  global.plot     = FALSE,
                  split.in        ="SH1")

## Every diagram its own scale ------------------------------
plots.by.state(outfile      = "plot.by.state.pdf",
               sp           = sp,
               seed_icu     = NULL,
               seed_dea     = NULL,
               iol          = iol,
               pspace       = pspace,
               rr           = rr,
               region       = "nuts2",
               fix.lim      = FALSE,
               filtered     = FALSE,
               Sec.Axis     = c("RMS"),
               silent       = FALSE,
               fk.cases     = rep(1/ 7,  7),
               fk.sec       = rep(1/15, 15),
               ind.states   = c(5),
               split.in     ="SH1")

################################################################################


setwd("../output/")
date.str<-"20211005"
icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")
dead       <- read.table(paste0(date.str,"dead_cases.csv"),head=TRUE,sep=",")
immune     <- read.table(paste0(date.str,"immune_cases.csv"),head=TRUE,sep=",")
icu$X<-1
rr<-list("healthy"=healthy,"inf_noncon"=inf_noncon,
         "inf_contag"=inf_contag,"ill_contag"=ill_contag,
         "ill_ICU"=icu,"dead"=dead,"immune"=immune)

setwd("../output-100p/")
date.str<-"20211001"
icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")
dead       <- read.table(paste0(date.str,"dead_cases.csv"),head=TRUE,sep=",")
immune     <- read.table(paste0(date.str,"immune_cases.csv"),head=TRUE,sep=",")
icu$X<-2
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
                  ind.states      = c(5),
                  global.plot     = FALSE,
                  split.in = "X")

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
               ind.states   = c(5),
               split.in = "X")

quit()

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



setwd("../output/")
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

plots.by.country (outfile         = "plot.by.country.pdf",
                  sp              = sp,
                  seed_icu        = iol$icu.cases.by.country,
                  seed_dea        = iol$dead.cases.by.country,
                  iol             = iol,
                  pspace          = pspace,
                  rr              = rr,
                  ind.states      = c(1,2,3,4,5,7,6),
                  global.plot     = TRUE)

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
quit()


plots.by.state(outfile      = "plot.by.state.pdf",
               sp           = sp,
               seed_icu     = NULL,
               seed_dea     = NULL,
               iol          = iol,
               pspace       = pspace,
               rr           = rr,
               region       = "state",
               fix.lim      = FALSE,
               filtered     = FALSE,
               ind.states   = c(5),
               split.in = "model")
rr<-list("healthy"=healthy,"inf_noncon"=inf_noncon,
          "inf_contag"=inf_contag,"ill_contag"=ill_contag,"ill_ICU"=rr.il_ICU,"dead"=dead,"immune"=immune ))

