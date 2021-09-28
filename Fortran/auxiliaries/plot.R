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

cur.dir <- getwd()

setwd("/lustre/nec/ws3/ws/hpcralf-BIB/cosmic-projection/Results-v12.0-2021-09-14_21:11:18/")                

iol<-readRDS("input-v12.0-2021-09-14_21:11:18.RDS")
pspace<-readRDS("pspace-v12.0-2021-09-14_21:11:18.RDS")
sp<-readRDS("static.params-v12.0-2021-09-14_21:11:18.RDS")

setwd(cur.dir)

args = commandArgs(trailingOnly=TRUE)

date.str  <- args[1]

icu        <- read.table(paste0(date.str,"ill_ICU_cases.csv"),head=TRUE,sep=",")
ill_contag <- read.table(paste0(date.str,"ill_contag_cases.csv"),head=TRUE,sep=",")
inf_contag <- read.table(paste0(date.str,"inf_contag_cases.csv"),head=TRUE,sep=",")
inf_noncon <- read.table(paste0(date.str,"inf_noncon_cases.csv"),head=TRUE,sep=",")
healthy    <- read.table(paste0(date.str,"healthy_cases.csv"),head=TRUE,sep=",")

sp$time_n <- 199

plots.by.country (outfile         = "test.pdf",
                  sp              = sp,
                  seed_icu        = iol$icu.cases.by.country,
                  seed_dea        = iol$dead.cases.by.country,
                  iol             = iol,
                  pspace          = pspace,
                  rr              = list("healthy"=healthy,"inf_noncon"=inf_noncon,
                                         "inf_contag"=inf_contag,"ill_contag"=ill_contag,
                                         "ill_ICU"=icu),
                  ind.states      = c(1,2,3,4,5),
                  global.plot     = TRUE,
                  split.in = "SH4")

