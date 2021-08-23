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
src.dir <- "/zhome/academic/HLRS/hlrs/hpcralf/CoSMic-states/cosmic/CoSMic/R"
for ( i in dir(src.dir,pattern="*.R$") ) { source(paste(src.dir,i,sep="/")) }

setwd("/zhome/academic/HLRS/hlrs/hpcralf/CoSMic-states/CoSMic_build/output/")                

iol<-readRDS("Results-v12.0-2021-06-29_09:18:35/input-v12.0-2021-06-29_09:18:35.RDS")
pspace<-readRDS("Results-v12.0-2021-06-29_09:18:35/pspace-v12.0-2021-06-29_09:18:35.RDS")
sp<-readRDS("Results-v12.0-2021-06-29_09:18:35/static.params-v12.0-2021-06-29_09:18:35.RDS")

icu<-read.table("20210629ill_ICU_cases.csv",head=TRUE,sep=",")

    plots.by.country (outfile         = "test.pdf",
                      sp              = sp,
                      seed_icu        = iol$icu.cases.by.country,
                      seed_dea        = iol$dead.cases.by.country,
                      iol             = iol,
                      pspace          = pspace,
                      rr              = list(1,2,3,4,"ill_ICU"=icu),
                      ind.states      = c(5),
                      global.plot     = TRUE)
