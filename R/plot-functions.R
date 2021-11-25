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
##==============================================================================
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
#' Plot timelines accross the complete country
#' 
#' The function plots timelines accross the complete country. Either fully
#' aggregated with global.plot = TRUE or aggregated once across the first
#' column of the latin hypercube, across each direct parameter with more than
#' one value and once across the parameter set of the first directv parameter.
#'
#' @import pracma
#' @import RColorBrewer
#' @import ggplot2
#' 
#' @export
plots.by.country <- function(outfile, sp, seed_icu, seed_dea,
                             iol, pspace,
                             rr,  ind.states=NULL, global.plot,
                             x.min=NULL,x.max=NULL, relative=FALSE,
                             silent=FALSE, split.in = NULL, y.max=NULL,
                             prog=NULL) {

    if (is.null(x.max) & is.null(sp$time_n)) {
        warning(paste("is.null(x.max) & is.null(sp$time_n) holds TRUE.",
                      "Both will be reset to max. value in R0change + 1"))
        sp$time_n <- max(unlist(sp$R0change))+1
    }
    
    if (is.null(x.min)) x.min <- sp$seed_date
    if (is.null(x.max)) x.max <- sp$seed_date + sp$time_n

    if (class(x.min) != "Date") x.min <- sp$seed_date + x.min
    if (class(x.max) != "Date") x.max <- sp$seed_date + x.max
       
    if (is.null(split.in) ) {
    ## --------------------------------------------------------------------------
    ## Select parameters from pspace to group by
    ## Extract parameter types from paramter list ---------------------
    pspace.types <- lapply(pspace,function(x){x$type})

    ## Extract parameters of type dist --------------------------------
    pspace.params <- lapply(pspace[which(pspace.types %in% c("dist","distv"))],
                            function(x){x$param})

    ## Extract standard deviations of type dist from pspace -----------
    pspace.sd <- lapply(pspace[which(pspace.types %in% c("dist","distv"))],
                        function(x){x$sd})
    
    ## Number of elements to vary in parameters of type dist ----------
    n.dist    <-  sum(unlist(pspace.sd)>0)
    
    ## 1.) The first of dist with sd > 1. Since all dist are sampled by LHS
    ##     there is no need to group by the other ones.
    if ( n.dist > 0 ) {
            gcols <- names(unlist(pspace.params[[1]])[which(unlist(pspace.sd) >0)])[1]
    } else {
        gcols <- NULL
    }
    
    ## 2.) All direct with length(param) > 1
    direct.cols <- which(unlist(lapply(
        pspace[which(pspace.types=="direct")],function(x){length(x$param)})) > 1)
    if ( length( direct.cols ) > 0 ) {
        if (length(gcols) == 0) {
            gcols <- names(direct.cols)
        } else {
            gcols <- c(gcols,names(direct.cols))
        }
    }
    
    ## 3.) All directv with not constant param
    if (length(pspace[which(pspace.types=="directv")]) == 1) {
        
        var.lens<-unlist(lapply(lapply(
            pspace[[which(pspace.types=="directv")]]$param,unique),length))
        
        if  (length(gcols) == 0) {
            gcols <- names(which(var.lens==max(var.lens))[1])
        } else {
            gcols <- c(gcols,
                       names(which(var.lens==max(var.lens))[1]))
        }
        
    } else if (length(pspace[which(pspace.types=="directv")]) > 1) {
        print("Aggregated plotting by parameter is only supported for one directv parameter")
        print("Skipping directv plots Sorry")
    }
    
    ## 4.) All directl with not constant param
    if (length(pspace[which(pspace.types=="directl")]) == 1) {
        
        lhc.dat <- do.call(rbind,lapply(pspace[[which(pspace.types=="directl")]]$param,unlist))
        
        var.lens<-unlist(lapply(apply(lhc.dat,2,unique),length))
        
        if  (length(gcols) == 0) {
            gcols <- colnames(lhc.dat)[which(var.lens==max(var.lens))[1]]
        } else {
            gcols <- c(gcols,
                       colnames(lhc.dat)[which(var.lens==max(var.lens))[1]])
        }
        
    } else if (length(pspace[which(pspace.types=="directl")]) > 1) {
        print("Aggregated plotting by parameter is only supported for one directl parameter")
        print("Skipping directl plots Sorry")
    }
    
        split.by <- gcols
        print(paste("Split by:",gcols))
    } else {
        split.by <- split.in
    }
    
    ## x-Axis ---------------------------------------------------
    time <- seq(from=sp$seed_date,to=sp$seed_date+sp$time_n-1,by="days")
    
    ## Vlines ----------------------------------------------------
    week.marks <- sp$seed_date+unlist(lapply(sp$R0change,function(x){x[1]}))-1

    if (is.null(ind.states)) {
        ind.states.to.plot <- names(rr)
    } else {

        if ( class(ind.states) == "character" ) {
            if ( all(ind.states %in% names(rr)) ) {
                ind.states.to.plot <- ind.states
            } else {
                warning(paste("Not all results requested to be plotted were found in rr","\n",
                              "Found in rr :",paste0(names(rr),collapse=","),"\n",
                              "Requested   :",paste0(ind.states,collapse=",")))
                ind.states.to.plot <- ind.states[ ind.states %in% names(rr) ]
            }
        } else if ( class(ind.states) == "numeric" ) {
            ind.states.to.plot <- names(rr)[ind.states]
            sp.states.map <- sp$sp.states[ind.states]
            names(sp.states.map) <- ind.states.to.plot
        } else {
            error("Please provide sp.states as vector of character strings or numeric values")
        }
        
    }
       
    ## Init list of plots ---------------------------------------
    plt <- list()
    
    for ( ii in  ind.states.to.plot ) {
        
        ## ------------------------------------------------------------------------------
        ## Plot States ------------------------------------------------------------------
        ## ------------------------------------------------------------------------------
        
        state.data<-as.data.frame(rr[[ii]])

        if (is.null(split.by)) {
            state.list <- list(state.data)
        } else {
            state.list <- split(state.data, state.data[,split.by])
        }

        if  ( (!is.null(prog)) & (length(state.list) != length(prog)) ) {
            prog <- NULL
            warning(paste0("length(prog) is not equal length(state.list) ",
                           "found by splitting current data.frame of state ",ii,"\n",
                           "prog is set to NULL."))
            prog <- NULL
        }
        
        if (relative & (ii == "ill_ICU")) {
            main.title <- paste(sp$country,paste0(ii," [%]"),sep=" - ")
            y.title <- paste(ii,"[%]")
        } else {
            main.title <- paste(sp$country,sp.states.map[ii],sep=" - ")        
            y.title <- paste(sp.states.map[ii])
        }

        plt[[ii]] <- ggplot() + ggtitle(main.title) +
            labs(x="Time [days]", y=y.title)
        
        fill.col <- rainbow(length(state.list))

        state.sum <- lapply(state.list,function(x){
            x[,c(paste0("X",1:sp$time_n),"iter")] %>%
                group_by(iter) %>% summarize_all(sum)})
        
        state.min  <- lapply(state.sum,function(x){
            apply(x[,paste0("X",1:sp$time_n)],2,min)})
        
        state.mean <- lapply(state.sum,function(x){
            apply(x[,paste0("X",1:sp$time_n)],2,mean)})
        
        state.max  <- lapply(state.sum,function(x){
            apply(x[,paste0("X",1:sp$time_n)],2,max)})
  
        if (global.plot) {
            state.min  <- apply(as.data.frame(state.min) ,1,min)
            state.mean <- apply(as.data.frame(state.mean),1,mean)
            state.max  <- apply(as.data.frame(state.max) ,1,max)

            if (relative & (ii == "ill_ICU")) {
                state.min  <- state.min  / sum(iol$icu.cap$sum) *100
                state.mean <- state.mean / sum(iol$icu.cap$sum) *100
                state.max  <- state.max  / sum(iol$icu.cap$sum) *100
            }
            
            fill.dt <- paste(round(unique(state.list[[1]][1,split.by]),3),collapse=" - ")
            
            df.gg <- data.frame(time=time,
                                dt.min = state.min, dt.mean=state.mean, dt.max=state.max,
                                data = fill.dt)
            
            plt[[ii]] <- plt[[ii]] +
                geom_ribbon(data=df.gg,aes(x=time, ymax=dt.max, ymin=dt.min), fill="red", alpha=.15) +
                geom_line(data=df.gg,aes(x=time,y = dt.max),  colour = "red",linetype = "dashed") +
                geom_line(data=df.gg,aes(x=time,y = dt.min),  colour = "red",linetype = "dashed") +
                geom_line(data=df.gg,aes(x=time,y = dt.mean), colour = "red", size=1)
                
        } else {
            
            for ( pg in seq(length(state.list)) ) {

                if ( is.null(prog) ) {
                    fill.dt <- paste(round(unique(state.list[[pg]][1,split.by]),3),collapse=" - ")
                } else {
                    fill.dt <- prog[pg]
                }
                
                if (relative & (ii == "ill_ICU")) {
                    state.min[[pg]]  <- state.min[[pg]]  / sum(iol$icu.cap$sum) *100
                    state.mean[[pg]] <- state.mean[[pg]] / sum(iol$icu.cap$sum) *100
                    state.max[[pg]]  <- state.max[[pg]]  / sum(iol$icu.cap$sum) *100
                }
                
                df.gg <- data.frame(time=time,
                                    dt.min = state.min[[pg]], dt.mean=state.mean[[pg]], dt.max=state.max[[pg]],
                                    data = fill.dt)
                 
                plt[[ii]] <- plt[[ii]] +
                    geom_ribbon(data=df.gg,aes(x=time, ymax=dt.max, ymin=dt.min, fill=data), alpha=.15) +
                    geom_line(data=df.gg,aes(x=time,y = dt.max , colour = data),linetype = "dashed") +
                    geom_line(data=df.gg,aes(x=time,y = dt.min , colour = data),linetype = "dashed") +
                    geom_line(data=df.gg,aes(x=time,y = dt.mean, colour = data), size=1)
                
            }
        }

        # R0change marks ---------------------------------------------
        plt[[ii]] <- plt[[ii]] + 
            geom_vline(xintercept = week.marks, linetype="dotted") +
            geom_text(data=data.frame(date=week.marks,event=seq(1,length(week.marks))),
                      mapping=aes(x=date, y=0, label=event),
                      size=4, angle=90, vjust=-0.4, hjust=0)
        
        # Observerd Ill, ICU -----------------------------------------
        if ( (ii == "ill_ICU") & !is.null(seed_icu) ) {

            if (relative) {
                observed <- seed_icu
                observed$cases <- observed$cases / sum(iol$icu.cap$sum) *100
            } else {
                observed <- seed_icu
            }
                          
            plt[[ii]] <- plt[[ii]] +
                geom_line(data=observed,
                          aes(x=date,y=cases,color="Observed")) +
                theme(legend.title = element_blank())

            m.cs <- "black"
            names(m.cs) <- "Observed"

            if (!global.plot) {
                if ("prog.RKI" %in% prog) {
                    m.cs <- c(m.cs, prog.RKI="red")
                }
                if ("constant" %in% prog) {
                    m.cs <- c(m.cs, constant="deepskyblue")
                }

                if (length(m.cs) < (length(state.list)+1)) {
                    m.cs <- c(m.cs,fill.col)
                    names(m.cs)[names(m.cs)==""] <- round(as.numeric(names(state.list)),3)
                }

                plt[[ii]] <- plt[[ii]] +
                    scale_color_manual(values = m.cs)+
                    scale_fill_manual (values = m.cs)
            } else {
                plt[[ii]] <- plt[[ii]] +
                    scale_color_manual(values = m.cs)
            }
        }

        # Observerd dead cases ---------------------------------------
        if ( (ii == "Dead") & !is.null(seed_dea) ) {

            plt[[ii]] <- plt[[ii]] +
                geom_line(data=seed_dea,
                          aes(x=date,y=cases,color="Observed"),size=1) +
                scale_color_manual(values = c('Observed' = 'black')) +
                theme(legend.title = element_blank())
        }
        
        if ( !is.null(split.by) ) {
            plt[[ii]] <- plt[[ii]] + labs(fill = split.by, colour=NULL)
        }

        plt[[ii]] <- plt[[ii]] + xlim(x.min,x.max)

        if ( ! is.null(y.max) ) {
            plt[[ii]] <- plt[[ii]] + ylim(0,y.max)
        }
        
    }

    if ( ! silent ) {

        ## Init PDF -------------------------------------------------
        pdf(outfile,width=16, height=9, onefile=TRUE)
        
        for ( ii in  ind.states.to.plot ) {
            print(plt[[ii]])
        }
        dev.off()
    } else {
        return(plt)
    }
}

################################################################################
#' Plot timelines accross each state
#' 
#' The function plots timelines accross each state, aggregated once across the
#' first column in the latin hypercube, ance across each direct parameter with
#' more than one value and once across the parameter set of the first directv
#' parameter.
#'
#' @import grid
#' @import gridExtra
#' @import pracma
#' @import RColorBrewer
#' @import ggplot2
#' 
#' @export
plots.by.state <- function(outfile, sp, seed_icu, seed_dea, iol,
                           pspace, rr, region, fix.lim,
                           filtered = FALSE, fk.cases=rep(1/7, 7),
                           Sec.Axis = "RMS", fk.sec=rep(1/15, 15),
                           sec.text = FALSE,
                           ind.states=NULL, silent=FALSE, relative=FALSE,
                           split.in = NULL, y.max=NULL, prog=NULL) {

        ## -------------------------------------------------------------------------
        ## Select parpameters from pspace to group by
        
    if ( is.null(split.in) ) {
        ## Extract parameter types from paramter list ---------------------
        pspace.types <- lapply(pspace,function(x){x$type})
        
        ## Extract parameters of type dist --------------------------------
        pspace.params <- lapply(pspace[which(pspace.types %in% c("dist","distv"))],
                                function(x){x$param})
        ## Extract standard deviations of type dist from pspace -----------
        pspace.sd <- lapply(pspace[which(pspace.types %in% c("dist","distv"))],
                            function(x){x$sd})
        
        ## Number of elements to vary in parameters of type dist ----------
        n.dist    <-  sum(unlist(pspace.sd)>0)
        
        ## 1.) The first of dist with sd > 1. Since all dist are sampled by LHS
        ##     there is no need to group by the other ones.
        if ( n.dist > 0 ) {
            gcols <- names(unlist(pspace.params[[1]])[which(unlist(pspace.sd) >0)])[1]
        } else {
            gcols <- NULL
        }
        
        ## 2.) All direct with length(param) > 1
        direct.cols <- which(unlist(lapply(
            pspace[which(pspace.types=="direct")],function(x){length(x$param)})) > 1)
        if ( length( direct.cols ) > 0 ) {
            if (length(gcols) == 0) {
                gcols <- names(direct.cols)
            } else {
                gcols <- c(gcols,names(direct.cols))
            }
        }
        
        ## 3.) All directv with not constant param
        if (length(pspace[which(pspace.types=="directv")]) == 1) {
            
            var.lens<-unlist(lapply(lapply(
                pspace[[which(pspace.types=="directv")]]$param,unique),length))
            
            if  (length(gcols) == 0) {
                gcols <- names(which(var.lens==max(var.lens))[1])
            } else {
                gcols <- c(gcols,
                           names(which(var.lens==max(var.lens))[1]))
            }
            
        } else if (length(pspace[which(pspace.types=="directv")]) > 1) {
            print("Aggregated plotting by parameter is only supported for one directv parameter")
            print("Skipping directv plots Sorry")
        }
        
        ## 4.) All directl with not constant param
        if (length(pspace[which(pspace.types=="directl")]) == 1) {
            
            lhc.dat <- do.call(rbind,lapply(pspace[[which(pspace.types=="directl")]]$param,unlist))
            
            var.lens<-unlist(lapply(apply(lhc.dat,2,unique),length))
            
            if  (length(gcols) == 0) {
                gcols <- colnames(lhc.dat)[which(var.lens==max(var.lens))[1]]
            } else {
                gcols <- c(gcols,
                           colnames(lhc.dat)[which(var.lens==max(var.lens))[1]])
            }
            
        } else if (length(pspace[which(pspace.types=="directl")]) > 1) {
            print("Aggregated plotting by parameter is only supported for one directl parameter")
            print("Skipping directl plots Sorry")
        }
        
        split.by <- gcols
        print(gcols)
    } else {
        split.by <- split.in
    }
    
    ## x-Axis ---------------------------------------------------
    time <- seq(from=as.Date(sp$seed_date),to=as.Date(sp$seed_date)+sp$time_n-1,by="days")

    ## Vlines ----------------------------------------------------
    week.marks <- as.Date(sp$seed_date)+unlist(lapply(sp$R0change,function(x){x[1]}))-1
     
    if (is.null(ind.states)) {
        ind.states.to.plot <- names(rr)
        sp.states.map <- ind.states.to.plot
        names(sp.states.map) <- ind.states.to.plot
    } else {

        if ( class(ind.states) == "character" ) {
            if ( all(ind.states %in% names(rr)) ) {
                ind.states.to.plot <- ind.states
            } else {
                warning(paste("Not all results requested to be plotted were found in rr","\n",
                              "Found in rr :",paste0(names(rr),collapse=","),"\n",
                              "Requested   :",paste0(ind.states,collapse=",")))
                ind.states.to.plot <- ind.states[ ind.states %in% names(rr) ]
            }
        } else if ( class(ind.states) == "numeric" ) {
            ind.states.to.plot <- names(rr)[ind.states]
            sp.states.map <- sp$sp.states[ind.states]
            names(sp.states.map) <- ind.states.to.plot
        } else {
            error("Please provide sp.states as vector of character strings or numeric values")
        }
       
    }

    if (relative & ("ill_ICU" %in% ind.states.to.plot)) {
        sp.states.map["ill_ICU"] <- paste("ill_ICU","[%]")
    } 
    

    ## Warnig if incompatible elements on secondary y-axis ---
    if ( ("RMS" %in% Sec.Axis) &
         ( ("R0effect.daily" %in% Sec.Axis) | ("R0effect.daily" %in% Sec.Axis) ) ) {
        warning(paste(paste(Sec.Axis,collapse=","),
                      " are selected to be plotted on secondary axis.\n",
                      "This results in incompatible scalings.\n",
                      "The R0effects  will only be observable qualitatively!"))
    }
    
    plt.pages <- list()
    
    for ( jj in ind.states.to.plot ) {

        ## The list of plots ---
        plt <- list()
               
        ## ---------------------------------------------------------------------
        ## Plot States ---------------------------------------------------------
        ## ---------------------------------------------------------------------
        if (region == "state") {
            region.ids <- unique(as.integer(unique(rr[[1]]$x.dist_id)/1000))
        } else if ( region == "nuts2" ) {
            region.ids <- unique(iol$counties[iol$counties$dist_id %in%
                                              unique(rr[[1]]$x.dist_id),"Nuts2"])
        } else {
            stop("Sorry region type ",region,"is not supported.!")
        }
        
        glob.max <- 0
        
        ii <- 1
        
        for ( st in region.ids ) { 
            
            if (region == "state") {
                ## Result data by state ----------------------------------------
                state.data <- as.data.frame(
                    rr[[jj]][as.integer(rr[[jj]]$x.dist_id/1000)==st,])

                ## init plot object --------------------------------------------
                plt[[ii]] <- ggplot() + ggtitle(paste(iol[["states"]][st,"Name"],
                                                  sp.states.map[jj],sep=" - "))
            }

            if ( region == "nuts2" ) {
                ## Result data by state ----------------------------------------
                state.data <- as.data.frame(
                    rr[[jj]][rr[[jj]]$x.dist_id %in%
                             iol[["counties"]][iol[["counties"]]$Nuts2 == st,"dist_id"],])
                ## init plot object --------------------------------------------
                plt[[ii]] <- ggplot() + ggtitle(
                                        paste(
                                            paste(unique(iol[["counties"]][iol$counties$Nuts2==st,
                                                                           c("Nuts2name","Nuts2")]),
                                                  collapse=" - "),
                                            sp.states.map[jj],sep=" - "))            
            }
                        
            if (is.null(split.by)) {
                state.list <- list(state.data)
            } else {
                state.list <- split(state.data, state.data[,split.by])
            }

            if ( (!is.null(prog)) & (length(state.list) != length(prog)) ) {
                prog <- NULL
                warning(paste0("length(prog) is not equal length(state.list) ",
                               "found by splitting current data.frame of state ",ii,"\n",
                               "prog is set to NULL."))
                prog <- NULL
            }
            
            plt[[ii]] <- plt[[ii]] + labs(x="Time [days]", y=sp.states.map[jj])
            
            fill.col <- rainbow(length(state.list))
            
            state.sum <- lapply(state.list,function(x) {
                x[,c(paste0("X",1:sp$time_n),"iter")] %>%
                    group_by(iter) %>% summarize_all(sum)})
            
            state.min  <- lapply(state.sum,function(x){
                apply(x[,paste0("X",1:sp$time_n)],2,min)})
            
            state.mean <- lapply(state.sum,function(x){
                apply(x[,paste0("X",1:sp$time_n)],2,mean)})
            
            state.max  <- lapply(state.sum,function(x){
                apply(x[,paste0("X",1:sp$time_n)],2,max)})

            if (region =="state") {
                cap <- iol$icu.cap$sum[st] / 100
            } else if (region =="nuts2") {
                cap <- iol$icu.cap.nuts2$sum[iol$icu.cap.nuts2$nuts2==st] / 100
            } else {
                warning(paste("Region type",region,"is not supported"))
            }
            
            for ( pg in seq(length(state.list)) ) {

                if (relative & (jj == "ill_ICU")) {
                    state.min[[pg]]  <- state.min[[pg]]  / cap
                    state.mean[[pg]] <- state.mean[[pg]] / cap
                    state.max[[pg]]  <- state.max[[pg]]  / cap
                }

                if ( is.null(prog) ) {
                    fill.dt <- paste(round(unique(state.list[[pg]][1,split.by]),3),collapse=" - ")
                } else {
                    fill.dt <- prog[pg]
                }
                
                if (filtered) {

                    df.gg <- data.frame(time=time,
                                        dt.min  = as.numeric(
                                            stats::filter(state.min[[pg]],fk.cases,sides=2)),
                                        dt.mean = as.numeric(
                                            stats::filter(state.mean[[pg]],fk.cases,sides=2)),
                                        dt.max  = as.numeric(
                                            stats::filter(state.max[[pg]],fk.cases,sides=2)),
                                        data = fill.dt) # paste0("Split - ",split.by,"=",fill.dt,"pg-",pg))
                } else {
                    df.gg <- data.frame(time=time,
                                        dt.min  = state.min[[pg]],
                                        dt.mean = state.mean[[pg]],
                                        dt.max  = state.max[[pg]],
                                        data = fill.dt) # paste0("Split - ",split.by,"=",fill.dt,"pg-",pg))
                }
                
                plt[[ii]] <- plt[[ii]] +
                    geom_ribbon(data=df.gg,aes(x=time, ymax=dt.max, ymin=dt.min, fill=data), alpha=.15) +
                    geom_line(data=df.gg,aes(x=time,y = dt.max,  color = data),linetype = "dashed") +
                    geom_line(data=df.gg,aes(x=time,y = dt.min,  color = data),linetype = "dashed") +
                    geom_line(data=df.gg,aes(x=time,y = dt.mean, color = data), size=1)
                
                ## Store global maximum ------------------------------
                glob.max[ii] <- max(c(glob.max[ii],max(df.gg[,2:4])))
                
            }
            
            ## Extract current x- & y-axis max -----------------------
            x.range <- layer_scales(plt[[ii]])$x$range$range[2]
            y.range <- layer_scales(plt[[ii]])$y$range$range[2]
            ## Scaling coefficient for secondary axis ----------------
            if ("RMS" %in% Sec.Axis) {
                s.coeff = y.range/50
                ax.scale <- 51
            } else {
                s.coeff = max(y.range)
                ax.scale <- 1
            } 

            ## R0change marks ----------------------------------------
            plt[[ii]] <- plt[[ii]] +
                geom_vline(xintercept = week.marks, linetype="dotted") +
                geom_text(data=data.frame(date=week.marks,event=seq(1,length(week.marks))),
                          mapping=aes(x=date, y=0, label=event),
                          size=3, angle=90, vjust=-0.1, hjust=0)
            
            ## Add observerd Ill, ICU and difference against observed ----------
            if ( (jj == "ill_ICU") & !is.null(seed_icu) ) {

                ## Get reference data ------------------------------------------
                if ( region == "state" ) {
                    ref.data <- seed_icu[seed_icu$state==iol$states[st,"Shortcut"],]
                }
                if ( region == "nuts2" ) {
                    ref.data <- as.data.frame(seed_icu[seed_icu$Nuts2==st,])
                }

                ref.data$cases.orig <- ref.data$cases
                    
                if (filtered) {
                    ref.data$cases <- as.numeric(
                        stats::filter(ref.data$cases,fk.cases,sides=2))
                }

                if (relative) {
                    ref.data$cases <- ref.data$cases / cap
                }
                                
                ## Add line for observed data-----------------------------------                
                plt[[ii]] <- plt[[ii]] +
                    geom_line(data=ref.data,
                              aes(x=date,y=cases,color='Observed')) +
                    theme(legend.title = element_blank())

                ## Standard colours are hard to differentiate ------------------
                ## Black for Observed
                m.cs <- "black"
                names(m.cs) <- "Observed"
                
                ## if limits are not globaly fixed -----------------------------
                if ( !fix.lim ) {

                    ## ---------------------------------------------------------
                    ## Plot R0effect as given in model parameters --------------
                    if ( "R0effect" %in% Sec.Axis) {

                        for ( pg in seq(length(state.list)) ) {

                            if ( is.null(prog) ) {
                                fill.dt <- paste(round(unique(state.list[[pg]][1,split.by]),3),collapse=" - ")
                            } else {
                                fill.dt <- prog[pg]
                            }
                            
                            if ( region == "state" ) {
                                gg.R0e <- data.frame(date = sp$seed_date +
                                                         unlist(lapply(sp$R0change,
                                                                       function(x){
                                                                           as.integer(x[1]+(x[2]-x[1])/2)})),
                                                     R0e = as.numeric(
                                                         state.list[[pg]][1,paste0(iol$states[st,"Shortcut"],
                                                         (1:length(sp$R0change)))]),
                                                     name= paste("R0effect",fill.dt))
                            }
                            if ( region == "nuts2" ) {
                                gg.R0e <- data.frame(date = sp$seed_date +
                                                         unlist(lapply(sp$R0change,
                                                                       function(x){
                                                                           as.integer(x[1]+(x[2]-x[1])/2)})),
                                                     R0e = as.numeric(
                                                         state.list[[pg]][1,paste0(st,
                                                         (1:length(sp$R0change)))]),
                                                     name= paste("R0effect",fill.dt))
                            }
                            
                            
                            ## Add line for R0effect on secondary y-axis -------
                            plt[[ii]] <- plt[[ii]] +
                                geom_line (data=gg.R0e,aes(x=date,y=R0e*!!enquo(s.coeff),color=name),
                                           linetype = "dashed",size=1) +
                                geom_point(data=gg.R0e,aes(x=date,y=R0e*!!enquo(s.coeff)))

                        }

                        tm.cs <- seq(length(state.list))
                        if (!is.null(prog)) {
                            names(tm.cs) <- paste(paste("R0effect",prog))
                        } else {
                            names(tm.cs) <- names(state.list)
                        }
                        
                        m.cs <- c(m.cs, tm.cs)

                    }
                    
                    ## ---------------------------------------------------------
                    ## Plot R0effect as used by model on daily basis -----------
                    if ( ("R0effect.daily" %in% Sec.Axis) & ( region == "state" ) ) {

                        for ( pg in seq(length(state.list)) ) {
                            ## Spread Weekly R0effect to get daily data ------------
                            R0effect.d  <- rep(as.numeric(
                                                     state.list[[pg]][1,paste0(iol$states[st,"Shortcut"],
                                                     (1:length(sp$R0change)))]),
                                               unlist(lapply(sp$R0change,function(x){x[2]-x[1]+1})))
                            
                            ## Smooth daily R0effects ------------------------------
                            R0effect.ds <- attenuate(R0effect.d,steps=sp$R0delay_days,type=sp$R0delay_type)
                            
                            gg.R0e.ds <- data.frame(time=time[1:length(R0effect.ds)],
                                                    R0e.ds=R0effect.ds,
                                                    name= paste0("Daily R0effect pg-",pg))
                            
                            ## Add line for daily R0effect on secondary y-axis -----
                            plt[[ii]] <- plt[[ii]] +
                                geom_line (data=gg.R0e.ds,aes(x=time,y=R0e.ds*!!enquo(s.coeff),color=name))
                        }
                        
                        tm.cs <- seq(length(state.list))
                        if (!is.null(prog)) {
                            names(tm.cs) <- paste(paste("R0effect.daily",prog))
                        } else {
                            names(tm.cs) <- names(state.list)
                        }
                        
                        m.cs <- c(m.cs, tm.cs)

                        
                    }
                    if ( ("R0effect.daily" %in% Sec.Axis) & ( region == "nuts2" ) ) {

                        warning("R0effect.daily for nuts2 not yet implemented")
                       
                    }                    

                    ## ---------------------------------------------------------
                    ## Plot first derivative of observed data ------------------
                    if ( "Observed.diff" %in% Sec.Axis) { 
                        
                        ## Smooth inner data with two sided kernel -------------
                        ref.data$cases <- as.numeric(
                            stats::filter(ref.data$cases.orig,fk.sec,sides=2))
                        
                        ## Smooth leftmost data with left truncated kernels ----
                        ref.data$cases[1:as.integer((length(fk.sec)-1)/2)] <- unlist(
                            lapply(
                                c((as.integer(length(fk.sec)/2)+1):(length(fk.sec)-1)),
                                function(x){sum(ref.data$cases.orig[1:x]*rep(1/x,x))}
                            ))
                        
                        ## Smooth rightmost data with rigth truncated kernels ----
                        len <- dim(ref.data)[1]
                        ref.data$cases[(dim(ref.data)[1]-ceiling((length(fk.sec)-1)/2)+1):
                                       dim(ref.data)[1]] <- unlist(
                            lapply(
                                c((length(fk.sec)-1):(as.integer((length(fk.sec)-1)/2)+1)),
                                function(x){sum(ref.data$cases.orig[(len+1-x):len]*rep(1/x,x))}
                            )) 
                        
                        ref.data$diff[2:dim(ref.data)[1]] <- c(ref.data$cases[2:dim(ref.data)[1]]-
                                                               ref.data$cases[1:dim(ref.data)[1]-1])
                        ref.data$diff[1] <- ref.data$diff[2] 

                        diff.weekly <- data.frame(date = sp$seed_date +
                                                      unlist(lapply(sp$R0change,
                                                                    function(x){as.integer(x[1]+(x[2]-x[1])/2)})),
                                                  diff = unlist(lapply(sp$R0change,function(x){
                                                      mean(ref.data[ref.data$date %in%
                                                                    seq(sp$seed_date+x[1],
                                                                        sp$seed_date+x[2],
                                                                        by="days"),"diff"])}
                                                      )))
                        
                        max.R0w <- max(pspace$R0effect$param$R0effect.ps1[,st])
                        max.dw <- max(diff.weekly$diff[!is.na(diff.weekly$diff)])
                        
                        diff.weekly$diff <- (diff.weekly$diff/max.dw+1)*max.R0w/2
                        diff.weekly$date <- diff.weekly$date - 13
                        
                        #rdd.max <- max(ref.data$diff)
                        #rdd.min <- min(ref.data$diff)#

                        #ref.data$diff <- ( ref.data$diff / rdd.max.val / 2 + 0.5 ) * s.coeff

                        ## Add line for gradient of observed data on y-axis) --------
                        plt[[ii]] <- plt[[ii]] +
                            geom_line(data=diff.weekly, #ref.data,
                                      aes(x=date,y=diff*!!enquo(s.coeff),
                                          color='Scaled.weekly.Grad.')) +
                            theme(legend.title = element_blank())
                                                tm.cs <- seq(length(state.list))

                        m.cs <- c(m.cs, Scaled.weekly.Grad.='Scaled weekly Grad.')

                        ## Standard colours are hard to differentiate ------------------
                        ## plt[[ii]] <- plt[[ii]] + scale_color_manual(
                        ##                              values = c("Observed"  = 'black',
                        ##                                         "mu"        = 'deepskyblue',
                        ##                                         "mu log."   = 'yellowgreen',
                        ##                                         "Scaled weekly Grad."  = 'darkorchid1'))
                    }
                    
                    ## ---------------------------------------------------------
                    ## Plot Root mean squares of observed against simulated ----
                    ## data and min against max of simulated data --------------
                    if ( "RMS" %in% Sec.Axis) {
                        
                        ## Generate deviation data against observed ------------
                        diff.i  <- data.frame(ref.data[ref.data$date %in% df.gg$time,c("date","cases")],
                                              diff=(df.gg[df.gg$time %in% ref.data$date,"dt.mean"] - 
                                                    ref.data[ref.data$date %in% df.gg$time,"cases"]  ),
                                              b.diff=(df.gg[df.gg$time %in% ref.data$date,"dt.max"] -
                                                      df.gg[df.gg$time %in% ref.data$date,"dt.min"] ))
                        
                        ## Aggregate deviation against observed to R0change timeframes
                        cmp <- data.frame(
                            do.call(
                                rbind,
                                ## Per R0change timeframe ------------------------------
                                lapply(sp$R0change,
                                       function(x) {
                                           ## 10% of mean of cases ---------------------
                                           cases <- mean(diff.i[diff.i$date %in% time[x[1]:x[2]],"cases"])*0.1
                                           ## Differences ------------------------------
                                           df    <- diff.i[diff.i$date %in% time[x[1]:x[2]],"diff"]
                                           ## Differences of solution bounds -----------
                                           bdf   <- diff.i[diff.i$date %in% time[x[1]:x[2]],"b.diff"]
                                           ## Root mean square against observed --------
                                           rms   <- sqrt(sum(df*df)/(x[2]-x[1]+1))
                                           ## Root mean square of bounds ---------------
                                           rmsb  <- sqrt(sum(bdf*bdf)/(x[2]-x[1]+1))
                                           
                                           df <- mean(df)
                                           
                                           round(c(cases,df,rms,rmsb),1)
                                       })),
                            ## Date values per R0change timeframe ----------------------
                            date = sp$seed_date + unlist(lapply(sp$R0change,
                                                                function(x){as.integer(x[1]+(x[2]-x[1])/2)})))
                        
                        names(cmp) <- c("Cases","Mean diff.","RMS","RMSB","date")
                        
                        ## Add line for 10% of cases per week (on second y-axis) ---
                        plt[[ii]] <- plt[[ii]] +
                            geom_line (data=cmp,aes(x=date,y=Cases*!!enquo(s.coeff),color='Scaled.cases.10pc'),
                                       linetype = "dashed", size = 1) +
                            geom_point(data=cmp,aes(x=date,y=Cases*!!enquo(s.coeff)))
                        if (sec.text) {
                            plt[[ii]] <- plt[[ii]] +
                                geom_text (data=cmp,aes(x=date,y=Cases*!!enquo(s.coeff),label=Cases),
                                           hjust=-0.1, vjust=-0.5)
                        }
                        
                        ## Add line for root mean square (on second y-axis) --------
                        plt[[ii]] <- plt[[ii]] +
                            geom_line (data=cmp,aes(x=date,y=RMS*!!enquo(s.coeff),color='RMS.cases')) +
                            geom_point(data=cmp,aes(x=date,y=RMS*!!enquo(s.coeff)))
                        if (sec.text) {
                            plt[[ii]] <- plt[[ii]] +
                                geom_text (data=cmp,aes(x=date,y=RMS*!!enquo(s.coeff),label=RMS),
                                           hjust=-0.1, vjust=-0.5)
                        }
                        
                        ## Add line for root mean square of bounds (on second y-axis) --------
                        plt[[ii]] <- plt[[ii]] +
                            geom_line (data=cmp,aes(x=date,y=RMSB*!!enquo(s.coeff),color='RMS.bounds')) +
                            geom_point(data=cmp,aes(x=date,y=RMSB*!!enquo(s.coeff)))
                        if (sec.text) {
                            plt[[ii]] <- plt[[ii]] +
                                geom_text (data=cmp,aes(x=date,y=RMSB*!!enquo(s.coeff),label=RMSB),
                                           hjust=-0.1, vjust=-0.5)
                        }

                        m.cs <- c(m.cs,
                                  Scaled.cases.10pc='Scaled.cases.10pc',
                                  RMS.cases='RMS.cases',
                                  RMS.bounds='RMS.bounds')
                    }

                    if ( ! is.null(Sec.Axis) ) {
                        ## Scale secondary y-axis ----------------------------------
                        plt[[ii]] <- plt[[ii]] +
                            scale_y_continuous(
                                sec.axis=sec_axis(~./max(.)*ax.scale,
                                                  name="-"))+
                            ## Zoom in the plot to original scale ------------------
                            coord_cartesian(xlim=c(sp$seed_date,as.POSIXct.Date(x.range)),
                                            ylim=c(0,y.range))
                    }

                    m.cs[-1] <- rainbow(length(m.cs)-1)
                } 
                
                if ("prog.RKI" %in% prog) {
                    m.cs <- c(m.cs, prog.RKI="red")
                }
                if ("constant" %in% prog) {
                    m.cs <- c(m.cs, constant="deepskyblue")
                }

                if (length(m.cs) < (length(state.list)+1)) {
                    m.cs <- c(m.cs,fill.col)
                    if (is.null(prog)) {
                        names(m.cs)[names(m.cs)==""] <- round(as.numeric(names(state.list)),3)
                    } else {
                        names(m.cs)[names(m.cs)==""] <- prog
                    }
                }
                
                plt[[ii]] <- plt[[ii]] +
                    scale_color_manual(values = m.cs)+
                    scale_fill_manual (values = m.cs)                  
            }

            ## Observerd dead cases ----------------------------------
            if ( (jj == "dead_cases") & !is.null(seed_dea) ) {
                if ( region == "state" ) {
                    plt[[ii]] <- plt[[ii]] +
                        geom_line(data=seed_dea[seed_dea$state==iol$states[st,"Name"],],
                                  aes(x=date,y=deaths,color='Observed'),size=1) +
                        scale_color_manual(values = c('Observed' = 'black')) +
                        theme(legend.title = element_blank())
                }
            }
            
            if ( !is.null(split.by) ) {
                plt[[ii]] <- plt[[ii]] + labs(fill = split.by, colour=NULL)
            }

            ii  <- ii + 1
            
            glob.max[ii] <- 0
        }
        
        if (fix.lim) {
            if ( ! is.null(y.max) ) { glob.max  <-  y.max }
            ii <- 1
            for ( st in region.ids ) {
                plt[[ii]] <- plt[[ii]] + ylim(0,max(glob.max))
                ii <- ii + 1
            }
        }
        plt.pages[[jj]] <- plt
    }

    if (!silent) {
        
        ## Init PDF -------------------------------------------------
        if (region == "state") {
            pdf(outfile,width=48, height=27, onefile=TRUE)
        } else if ( region == "nuts2" ) {
            pdf(outfile,width=64, height=36, onefile=TRUE)
        } else {
            stop("Sorry region type ",region,"is not supported.!")
        }
        
        for ( jj in ind.states.to.plot ) {
            grid.arrange(grobs=plt.pages[[jj]], ncol=round(sqrt(length(region.ids))), legend="TRUE",
                         top = textGrob(paste("Sample Size :",pspace[["sam_size"]]$param,
                                              "| Sim. Period :"  ,sp$seed_date,"-",
                                              sp$seed_date+sp$time_n,
                                              "| R0 :",pspace[["R0"]],
                                              "| inf_dur :",sp$inf_dur,
                                              "| cont_dur :",sp$cont_dur,
                                              "| ill_dur :",sp$ill_dur,
                                              "| icu_dur :",pspace[["icu_dur"]]$param,
                                              "| icu_per_day :",
                                              paste(sp$icu_per_day,sep="",collapse=","),
                                              "| w_int :",pspace[["w_int"]]$param
                                              ),
                                        gp=gpar(fontsize=20,font=3)))  
        }
        
        dev.off()
    } else {
        return(plt.pages[[jj]])
    }
    
}


################################################################################
#' Plot R0effects over R0changes    
#' 
#' The function plots timelines of the R0effects per state or NUTS2 region.
#'
#' @export
plot.R0effect <- function(R0effect,sp,outfile=NULL,silent=FALSE) {

    dt   <-  sp$seed_date + unlist(lapply(sp$R0change,function(x){as.integer(x[1]+(x[2]-x[1])/2)}))
    d2 <- dim(R0effect)[2]
    if (is.null(d2)) d2 <- 1
    
    date <- rep(dt,d2)

    if (d2 > 1) {
        gg <- data.frame(date=rep(date,d2), R0effect=unlist(R0effect[1:length(sp$R0change),]),
                         region=rep(names(R0effect),each=length(sp$R0change)))
    } else {
        gg <- data.frame(date=date[1:length(sp$R0change)],
                         R0effect=R0effect[1:length(sp$R0change)],
                         region=rep("R0effect",each=length(sp$R0change)))
    }
    labs <- data.frame(date=dt,dat=rep(1.0,length(dt)),lab=seq(1,length(dt)))
    
    plt <- ggplot(data=gg) +
        geom_line (aes(x=date,y=R0effect,color=region)) +
        geom_point(aes(x=date,y=R0effect,color=region)) +
        geom_point(data=labs,
                   aes(x=date,y=dat)) +
        geom_text (data=labs,
                   aes(x=date,y=dat,label=lab),
                   hjust=0.5, vjust=-1.0) +
        ylim(0,1)

    if(!silent) {
        if (is.null(outfile)) {
            print(plt)
        } else {
            pdf(outfile,width=16, height=9)
            print(plt)
            dev.off()
        }
    }

    return(plt)
}

plot.states.combined <- function(outfile, sp, seed_icu, seed_dea, iol,
                                 pspace, rr, region, ind.states) {
       
    plt <- list()

    ## ---------------------------------------------------------------------
    ## Plot States ---------------------------------------------------------
    ## ---------------------------------------------------------------------
    if (region == "state") {
        region.ids <- unique(as.integer(unique(rr[[1]]$x.dist_id)/1000))
    } else if ( region == "nuts2" ) {
        region.ids <- unique(iol$counties[iol$counties$dist_id %in%
                                          unique(rr[[1]]$x.dist_id),"Nuts2"])
    } else {
        stop("Sorry region type ",region,"is not supported.!")
    }

    sp.states  <- c("Healthy / Never infected",    
                    "Infected, not contagious",        
                    "Infected, contagious",            
                    "Ill, contagious",                 
                    "Ill, ICU",                        
                    "Immune",                          
                    "Dead")
        
    ## Plot non filtered --------------------------------------------
    plt.nf <- plots.by.state(paste0(outfile,"nf.pdf"), sp, seed_icu, seed_dea, iol,
                             pspace, rr, region,
                             fix.lim  = FALSE,
                             filtered = FALSE,
                             ind.states=ind.states,
                             silent=TRUE)

    ## Plot filtered ------------------------------------------------
    plt.f  <- plots.by.state(paste0(outfile,"f.pdf"), sp, seed_icu, seed_dea, iol,
                             pspace, rr, region,
                             fix.lim  = FALSE,
                             filtered = TRUE,
                             ind.states=ind.states,
                             silent=TRUE)

    ## Plot filtered with R0 ----------------------------------------
    plt.R  <- plots.by.state(paste0(outfile,"f1.pdf"), sp, seed_icu, seed_dea, iol,
                             pspace, rr, region,
                             fix.lim  = FALSE,
                             filtered = TRUE,
                             Sec.Axis=c("R0effect","R0effect.daily","Observed.diff"),
                             ind.states=ind.states,
                             silent=FALSE,
                             fk.cases=rep(1/7, 7),
                             fk.sec=rep(1/15, 15))
    ii <- 1
    plt <- list()
    for ( st in region.ids ) {
        plt[[ii]] <- grid.arrange(
            grobs=list(plt.nf[[ii]],plt.f[[ii]],plt.R[[ii]]),
                       ncol=1)
        ii <- ii+1
    }

    ## Init PDF -------------------------------------------------
    if (region == "state") {
        pdf(paste0(outfile,".pdf"),width=48, height=81, onefile=TRUE)
    } else if ( region == "nuts2" ) {
        pdf(paste0(outfile,".pdf"),width=64, height=96, onefile=TRUE)
    } else {
        stop("Sorry region type ",region,"is not supported.!")
    }
    
    grid.arrange(grobs=plt,ncol=round(sqrt(length(region.ids))))

    dev.off()
    
}

