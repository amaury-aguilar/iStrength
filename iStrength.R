#!/usr/bin/env Rscript
# iStrength.v11.R
# Calculo de enriquecimientos...
# Author: Oscar Amaury Aguilar Lomas
# Date: May 5, 2025
# Uso:
# Rscript iStrength.v9.beta.R archivo_de_configuracion nombre_general
# Corregido problema con el ID de los archivos
#Bibliotecas necesarias -----------
#install.packages("data.table")
library(data.table)
library(parallel)
library(matrixStats)
#Si no está instalado, este es el comando para instalar straw :
#remotes::install_github("aidenlab/straw/R")
#- Variables ----------------------
variables <- commandArgs(trailingOnly = TRUE)
configuracion <- read.table(as.character(variables[1]),header = T,row.names = 1)
#configuracion <- read.table("config/loops_TADs.3925MvP.iStrength.config",
 #                           header = T,row.names = 1)
#configuracion <- read.table("config/pyHICCUPS.3kb_4kb_5kb_6kb_8kb_10kb.iStrength.config",header = T,row.names = 1)
colnames(configuracion[,1:ncol(configuracion)-1])
# i <- "TADs_3925MvP_Ctcf_WT_test"
for (i in colnames(configuracion[,1:ncol(configuracion)-1])){
  nomen <- i
  cores <- as.character(configuracion["cores:", 1])
  chr.list.var <- as.character(configuracion["genome_reference:", 1])
  (hic.dir <- as.character(configuracion["hic_directories:", 1]))
  (hic.file <- as.character(configuracion["hic_file:", i]))
  (regiones_dir <- as.character(configuracion["regiones_dir:", 1]))
  (peak.file.var <- as.character(configuracion["regiones:", i]))
  (resol <- as.numeric(as.character(configuracion["resolucion:", 1])))
  (choose_window <- as.character(configuracion["choose_window:", 1]))
  (window <- as.numeric(as.character(configuracion["window:", 1])))
  (window_ratio <- as.numeric(as.character(configuracion["window_ratio:", 1])))
  (dmin <- as.numeric(as.character(configuracion["dmin:", 1])))
  (window1 <- as.numeric(as.character(configuracion["window_anchors_in:", 1])))
  (window2 <- as.numeric(as.character(configuracion["window_anchors_out:", 1])))
  (stripe_w <- as.numeric(as.character(configuracion["stripe_width:", 1])))
  (id <- as.character(configuracion["id:", 1]))
  # Windows - Load-chr-table----------------------------------------
  (window <- floor(window/resol)*resol)
  (stripe_width <- floor(stripe_w/resol))
  #chr.list <- read.table(paste0("genomes/",chr.list.var))
  # Load-chr-table----------------------------------------
  mm10 <- data.table(chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"),
                     size = c(195471971, 182113224, 160039680, 156508116, 151834684,149736546, 145441459, 129401213, 124595110, 130694993,122082543, 120129022, 120421639, 124902244, 104043685,98207768, 94987271, 90702639, 61431566, 171031299, 91744698))
  mm10_test <- data.table(chr = c("18", "19"),size = c(90702639, 61431566))
  hg38 <- data.table(chr=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"),
                     size=c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415))
  
  genome.list <- list(mm10=mm10,hg38=hg38,mm10_test=mm10_test)
  
  chr.list <- genome.list[[chr.list.var]]
  
  #Take time and Set dir -------------------------------------------------------
  start.time <- Sys.time()
  options(scipen=999)
  dir <- getwd()
  setwd(dir)
  #Take time and Set dir -------------------------------------------------------
  dir.create(file.path(dir,"iStrength"))
  dir.create(file.path(dir,"iStrength", nomen))
  dir.create(file.path(dir,"iStrength", nomen, resol))
  # Meta #------------------------------------------------------------
  domains.bed <- data.table(chr=character(),
                            start=numeric(),end=numeric(),
                            id=numeric(),
                            length=numeric(),
                            outTAD_win=numeric(),
                            upTAD_start_obs=numeric(),
                            downTAD_start_obs=numeric(),
                            upTAD_end_obs=numeric(),
                            downTAD_end_obs=numeric(),
                            interTAD_start_obs=numeric(),
                            iStrength_start_obs=numeric(),
                            interTAD_end_obs=numeric(),
                            iStrength_end_obs=numeric(),
                            interTAD_start_log=numeric(),
                            iStrength_start_log=numeric(),
                            interTAD_end_log=numeric(),
                            iStrength_end_log=numeric(),
                            interTAD_start_median=numeric(),
                            iStrength_start_median=numeric(),
                            interTAD_end_median=numeric(),
                            iStrength_end_median=numeric(),
                            upTAD_start_oe=numeric(),
                            downTAD_start_oe=numeric(),
                            upTAD_end_oe=numeric(),
                            downTAD_end_oe=numeric(),
                            interTAD_start_oe=numeric(),
                            iStrength_start_oe=numeric(),
                            interTAD_end_oe=numeric(),
                            iStrength_end_oe=numeric(),
                            interTAD_start_oe_median=numeric(),
                            iStrength_start_oe_median=numeric(),
                            interTAD_end_oe_median=numeric(),
                            iStrength_end_oe_median=numeric(),
                            tad_counts=numeric(),
                            tad_obs=numeric(),
                            tad_o_avg=numeric(),
                            tad_log=numeric(),
                            tad_log_avg=numeric(),
                            tad_oe=numeric(),
                            loop_counts=numeric(),
                            loop_obs=numeric(),
                            loop_o_avg=numeric(),
                            loop_log=numeric(),
                            loop_log_avg=numeric(),
                            loop_oe=numeric(),
                            loop_ratio_obs=numeric(),
                            loop_ratio_oe=numeric(),
                            up_stripe_counts=numeric(),
                            up_stripe_obs=numeric(),
                            up_stripe_o_avg=numeric(),
                            up_stripe_log=numeric(),
                            up_stripe_log_avg=numeric(),
                            up_stripe_oe=numeric(),
                            down_stripe_counts=numeric(),
                            down_stripe_obs=numeric(),
                            down_stripe_o_avg=numeric(),
                            down_stripe_log=numeric(),
                            down_stripe_log_avg=numeric(),
                            down_stripe_oe=numeric(),
                            up_stripe_ratio_counts=numeric(),
                            up_stripe_ratio_obs=numeric(),
                            up_stripe_ratio_oe=numeric(),
                            down_stripe_ratio_counts=numeric(),
                            down_stripe_ratio_obs=numeric(),
                            down_stripe_ratio_oe=numeric(),
                            avg=numeric())
  # Functions ------------------------------------------------------------------
  matistics <- function(i){
    #Start boundary ------------------------------------------------------------
    outTAD_win <-tads.chr[i,6]
    start_domain <- tads.chr[i,2]
    end_domain <- tads.chr[i,3]
    
    upTAD.start.range <- seq.int(start_domain-outTAD_win,start_domain,resol)
    upTAD.start.intervalo <- 
      data.table(x = rep(upTAD.start.range,
                         times = length(upTAD.start.range)),
                 y = rep(upTAD.start.range,
                         each = length(upTAD.start.range)))[x < y]
    setkey(upTAD.start.intervalo, x, y)
    upTAD.start.region <- mtz[upTAD.start.intervalo, on = .(x, y)]
    upTAD.start.region[is.na(upTAD.start.region)] <- 0
    upTAD.start <- colMeans(upTAD.start.region[,-c("x","y")])
    upTAD.start.median <- 
      colMedians(as.matrix(upTAD.start.region[,-c("x","y")]))
    
    downTAD.start.range <- seq.int(start_domain,start_domain+outTAD_win,resol)
    downTAD.start.intervalo <-
      data.table(x = rep(downTAD.start.range, 
                         times = length(downTAD.start.range)),
                 y = rep(downTAD.start.range, 
                         each = length(downTAD.start.range)))[x < y]
    setkey(downTAD.start.intervalo, x, y)
    downTAD.start.region <- mtz[downTAD.start.intervalo, on = .(x, y)]
    downTAD.start.region[is.na(downTAD.start.region)] <- 0
    downTAD.start <- colMeans(downTAD.start.region[,-c("x","y")])
    downTAD.start.median <- 
      colMedians(as.matrix(downTAD.start.region[,-c("x","y")]))
    
    interTAD.start.intervalo <- data.table(
      x = rep(head(upTAD.start.range,-1), 
              times = length(head(upTAD.start.range,-1))),
      y = rep(downTAD.start.range[-1], each = length(downTAD.start.range[-1]))
    )[x < y]
    setkey(interTAD.start.intervalo, x, y)
    length(interTAD.start.intervalo$x)
    interTAD.start.region <- mtz[interTAD.start.intervalo, on = .(x, y)]
    interTAD.start.region[is.na(interTAD.start.region)] <- 0
    interTAD.start <- colMeans(interTAD.start.region[,-c("x","y")])
    interTAD.start.median <- 
      colMedians(as.matrix(interTAD.start.region[,-c("x","y")]))
    
    #End boundary  -------------------------------------------------------------  
    upTAD.end.range <- seq.int(end_domain-outTAD_win,end_domain,resol)
    upTAD.end.intervalo <- data.table(x = rep(upTAD.end.range, 
                                              times = length(upTAD.end.range)),
                                      y = rep(upTAD.end.range,
                                              each = length(upTAD.end.range)))[x < y]
    setkey(upTAD.end.intervalo, x, y)
    upTAD.end.region <- mtz[upTAD.end.intervalo, on = .(x, y)]
    upTAD.end.region[is.na(upTAD.end.region)] <- 0
    upTAD.end <- colMeans(upTAD.end.region[,-c("x","y")])
    upTAD.end.median <- 
      colMedians(as.matrix(upTAD.end.region[,-c("x","y")]))
    
    downTAD.end.range <- seq.int(end_domain,end_domain+outTAD_win,resol)
    downTAD.end.intervalo <- data.table(x = rep(downTAD.end.range, 
                                                times = length(downTAD.end.range)),
                                        y = rep(downTAD.end.range,
                                                each = length(downTAD.end.range)))[x < y]
    setkey(downTAD.end.intervalo, x, y)
    downTAD.end.region <- mtz[downTAD.end.intervalo, on = .(x, y)]
    downTAD.end.region[is.na(downTAD.end.region)] <- 0
    downTAD.end <- colMeans(downTAD.end.region[,-c("x","y")])
    downTAD.end.median <- 
      colMedians(as.matrix(downTAD.end.region[,-c("x","y")]))
    
    interTAD.end.intervalo <- data.table(
      x = rep(head(upTAD.end.range,-1), times = length(head(upTAD.end.range,-1))),
      y = rep(downTAD.end.range[-1], each = length(downTAD.end.range[-1]))
    )[x < y]
    setkey(interTAD.end.intervalo, x, y)
    interTAD.end.region <- mtz[interTAD.end.intervalo, on = .(x, y)]
    interTAD.end.region[is.na(interTAD.end.region)] <- 0
    interTAD.end <- colMeans(interTAD.end.region[,-c("x","y")])
    interTAD.end.median <- 
      colMedians(as.matrix(interTAD.end.region[,-c("x","y")]))
    
    #iStrength ------------------------------------------------------  
    intraTAD.start <- upTAD.start+downTAD.start
    iStrength.start <- (intraTAD.start-interTAD.start)/(intraTAD.start+interTAD.start)
    intraTAD.start.median <- upTAD.start.median+downTAD.start.median
    iStrength.start.median <-
      (intraTAD.start.median-interTAD.start.median)/(intraTAD.start.median+interTAD.start.median)
    
    
    intraTAD.end <- upTAD.end+downTAD.end
    iStrength.end <- (intraTAD.end-interTAD.end)/(intraTAD.end+interTAD.end)
    intraTAD.end.median <- upTAD.end.median+downTAD.end.median
    iStrength.end.median <-
      (intraTAD.end.median-interTAD.end.median)/(intraTAD.end.median+interTAD.end.median)
    
    #TAD-------------------------------------------------------------------------
    tad.range <- seq.int(start_domain,end_domain,resol)
    tad.intervalo <- data.table(x = rep(tad.range, 
                                        times = length(tad.range)),
                                y = rep(tad.range,
                                        each = length(tad.range)))[x < y]
    setkey(tad.intervalo, x, y)
    tad.region <- mtz[tad.intervalo, on = .(x, y)]
    tad.region[is.na(tad.region)] <- 0
    tad.means <- colMeans(tad.region[,-c("x","y")])
    
    #Stripes-------------------------------------------------------------------------
    up_stripes.range <- 
      seq.int(start_domain,end_domain-(stripe_width*resol),resol)
    up_stripes.intervalo <- data.table(x = rep(up_stripes.range[1:stripe_width],
                                               each=length(up_stripes.range)),
                                       y = rep(up_stripes.range,stripe_width))[x < y]
    setkey(up_stripes.intervalo, x, y)
    up_stripes.region <- mtz[up_stripes.intervalo, on = .(x, y)]
    up_stripes.region[is.na(up_stripes.region)] <- 0
    up_stripe <- colMeans(up_stripes.region[,-c("x","y")])
    
    down_stripes.range <- seq.int(start_domain+(stripe_width*resol),end_domain,resol)
    down_stripes.intervalo <- data.table(x = rep(down_stripes.range,stripe_width),
                                         y = rep(tail(down_stripes.range,stripe_width),
                                                 each=length(down_stripes.range)))[x < y]
    setkey(down_stripes.intervalo, x, y)
    down_stripes.region <- mtz[down_stripes.intervalo, on = .(x, y)]
    down_stripes.region[is.na(down_stripes.region)] <- 0
    down_stripe <- colMeans(down_stripes.region[,-c("x","y")])
    # both stripes _   _   _   _   _   _   _   _   _   _
    anchor_up_stripes.range <- seq.int(start_domain,start_domain+(stripe_width*resol),resol)
    anchor_down_stripes.range <- seq.int(end_domain-(stripe_width*resol),end_domain,resol)
    anchor_stripes.intervalo <- data.table(x = rep(anchor_up_stripes.range,length(anchor_up_stripes.range)),
                                           y = rep(anchor_down_stripes.range,
                                                   each=length(anchor_down_stripes.range)))[x < y]
    anchor_stripes.region <- mtz[anchor_stripes.intervalo, on = .(x, y)]
    anchor_stripes.region[is.na(anchor_stripes.region)] <- 0
    stripes.region <- rbind(up_stripes.region,down_stripes.region)
    #TAD_NoStripes_NoLoop & Stripes Ratio ---------------------------
    tad.no_stripes.region <- tad.region[!stripes.region, on = .(x, y)]
    tad.no_stripes.no_loop.region <- tad.no_stripes.region[!anchor_stripes.region, on = .(x, y)]
    tad.no_stripes.no_loop.region[is.na(tad.no_stripes.no_loop.region)] <- 0
    tad.no_stripes.no_loop <- colMeans(tad.no_stripes.no_loop.region[,-c("x","y")])
    
    up_stripe.ratio <- up_stripe/tad.no_stripes.no_loop
    down_stripe.ratio <- down_stripe/tad.no_stripes.no_loop
    
    #Loop -------------------------------------------------------------------------
    loop_anchor_down <- seq.int(start_domain-window2,start_domain+window1,resol)
    loop_anchor_up <- seq.int(end_domain-window1,end_domain+window2,resol)
    loop.intervalo <- data.table(x = rep(loop_anchor_down, 
                                         times = length(loop_anchor_down)),
                                 y = rep(loop_anchor_up,
                                         each = length(loop_anchor_up)))[x < y]
    setkey(loop.intervalo, x, y)
    loop.region <- mtz[loop.intervalo, on = .(x, y)]
    loop.region[is.na(loop.region)] <- 0
    loop <- colMeans(loop.region[,-c("x","y")])
    loop.ratio <- loop/tad.no_stripes.no_loop
    
    #Data Frame ------------------------------------------------------------------
    domains.df <- data.table(chr=tads.chr[i,1],
                             start=start_domain,end=end_domain,
                             id=tads.chr[i,4],
                             length=tads.chr[i,5],
                             outTAD_win=tads.chr[i,6],
                             upTAD_start_obs=upTAD.start["obs"],
                             downTAD_start_obs=downTAD.start["obs"],
                             upTAD_end_obs=upTAD.end["obs"],
                             downTAD_end_obs=downTAD.end["obs"],
                             interTAD_start_obs=interTAD.start["obs"],
                             iStrength_start_obs=iStrength.start["obs"],
                             interTAD_end_obs=interTAD.end["obs"],
                             iStrength_end_obs=iStrength.end["obs"],
                             interTAD_start_log=interTAD.start["log"],
                             iStrength_start_log=iStrength.start["log"],
                             interTAD_end_log=interTAD.end["log"],
                             iStrength_end_log=iStrength.end["log"],
                             interTAD_start_median=interTAD.start.median["log"],
                             iStrength_start_median=iStrength.start.median["log"],
                             interTAD_end_median=interTAD.end.median["log"],
                             iStrength_end_median=iStrength.end.median["log"],
                             upTAD_start_oe=upTAD.start["oe"],
                             downTAD_start_oe=downTAD.start["oe"],
                             upTAD_end_oe=upTAD.end["oe"],
                             downTAD_end_oe=downTAD.end["oe"],
                             interTAD_start_oe=interTAD.start["oe"],
                             iStrength_start_oe=iStrength.start["oe"],
                             interTAD_end_oe=interTAD.end["oe"],
                             iStrength_end_oe=iStrength.end["oe"],
                             interTAD_start_oe_median=interTAD.start.median["oe"],
                             iStrength_start_oe_median=iStrength.start.median["oe"],
                             interTAD_end_oe_median=interTAD.end.median["oe"],
                             iStrength_end_oe_median=iStrength.end.median["oe"],
                             tad_counts=tad.means["counts"],
                             tad_obs=tad.means["obs"],
                             tad_o_avg=tad.means["o_avg"],
                             tad_log=tad.means["log"],
                             tad_log_avg=tad.means["log_avg"],
                             tad_oe=tad.means["oe"],
                             loop_counts=loop["counts"],
                             loop_obs=loop["obs"],
                             loop_o_avg=loop["o_avg"],
                             loop_log=loop["log"],
                             loop_log_avg=loop["log_avg"],
                             loop_oe=loop["oe"],
                             loop_ratio_obs=loop.ratio["obs"],
                             loop_ratio_oe=loop.ratio["oe"],
                             up_stripe_counts=up_stripe["counts"],
                             up_stripe_obs=up_stripe["obs"],
                             up_stripe_o_avg=up_stripe["o_avg"],
                             up_stripe_log=up_stripe["log"],
                             up_stripe_log_avg=up_stripe["log_avg"],
                             up_stripe_oe=up_stripe["oe"],
                             down_stripe_counts=down_stripe["counts"],
                             down_stripe_obs=down_stripe["obs"],
                             down_stripe_o_avg=down_stripe["o_avg"],
                             down_stripe_log=down_stripe["log"],
                             down_stripe_log_avg=down_stripe["log_avg"],
                             down_stripe_oe=down_stripe["oe"],
                             up_stripe_ratio_counts=up_stripe.ratio["counts"],
                             up_stripe_ratio_obs=up_stripe.ratio["obs"],
                             up_stripe_ratio_oe=up_stripe.ratio["oe"],
                             down_stripe_ratio_counts=down_stripe.ratio["counts"],
                             down_stripe_ratio_obs=down_stripe.ratio["obs"],
                             down_stripe_ratio_oe=down_stripe.ratio["oe"],
                             avg=avg
    )
    print(domains.df)
  }
  # Chr-Loop---------------------------------------------------------------
  # j <- "18"
  for (j in row.names(chr.list)) {
    (chr <- as.character(chr.list[as.integer(j),chr]))
    # Si hay loops ... ---------------------------- #
    print(paste0("Checking loops file for chr",chr))
    peaks <- read.table(paste0(regiones_dir,peak.file.var))
    peaks$V1 <- gsub("chr","",peaks$V1)
    peaks.chr.pre <- peaks[which(peaks$V1==chr),]
    if(length(peaks.chr.pre$V1)>0){
      
    bins.chr <- (ceiling(chr.list[as.integer(j),size]/resol)**2)/2
    #__Calculando____________________________
    print(paste0("Calculando fuerza de interaccion en anclas de loops para el chr",chr))
    #Cargando matrices ---------------------------- #
    print(paste0("Cargando matrices .hic de chr",chr))
    mtz <- strawr::straw("KR", paste0(hic.dir,hic.file), chr, chr, "BP", resol,matrix="observed")
    colnames(mtz) <- c("x","y","obs")
    mtz$log <- log10(mtz$obs+1)
    avg <- (sum(mtz$obs,na.rm=T))/bins.chr
    mtz$o_avg <- mtz$obs/avg
    mtz$log_avg <- log10(mtz$o_avg+1)
    mtz.oe <- strawr::straw("KR", paste0(hic.dir,hic.file), chr, chr, "BP", resol,matrix="oe")
    colnames(mtz.oe) <- c("x","y","oe")
    mtz.raw <- strawr::straw("NONE", paste0(hic.dir,hic.file), chr, chr, "BP", resol,matrix="observed")
    colnames(mtz.raw) <- c("x","y","counts")
    
    # Set_Keys
    setDT(mtz)
    setDT(mtz.raw)
    setDT(mtz.oe)
    setkey(mtz, x, y)
    setkey(mtz.raw, x, y)
    setkey(mtz.oe, x, y)
    mtz <- mtz.oe[mtz, on = .(x, y)]
    mtz <- mtz.raw[mtz, on = .(x, y)]
    
      # Preparando loops --------------
      if(id == "id"){
        tads.chr <- data.frame(chr=peaks.chr.pre$V1,
                               start=floor((peaks.chr.pre$V2)/as.numeric(resol))*resol,
                               end=floor((peaks.chr.pre$V3)/as.numeric(resol))*resol,
                               id=peaks.chr.pre$V4,
                               length=peaks.chr.pre$V3-peaks.chr.pre$V2)
      }else{
        tads.chr <- data.frame(chr=peaks.chr.pre$V1,
                               start=floor((peaks.chr.pre$V2)/as.numeric(resol))*resol,
                               end=floor((peaks.chr.pre$V3)/as.numeric(resol))*resol,
                               id=paste0(nomen,"_chr",chr,"_",seq(1,length(peaks.chr.pre$V1))),
                               length=peaks.chr.pre$V3-peaks.chr.pre$V2)
      }
      tads.chr <- tads.chr[tads.chr$length >= dmin,]
      
      if(choose_window=="ratio"){
        tads.chr$outTAD_win <- tads.chr$length*window_ratio
        window_label <- paste0(window_ratio,"_TADLength")
      }else if(choose_window=="fixed"){
        tads.chr$outTAD_win <- window
        window_label <- paste0(window/1000,"_kb")
        }
      
      conteos <- function(){
        print(paste0("Calculando fuerza de interacciones en dominio de ",chr))
        tad_parallel.start_time <- Sys.time()
        tads.statistics <- mclapply(1:length(tads.chr$id), matistics, mc.cores = cores)
        tads.statistics.dt <- rbindlist(lapply(tads.statistics, function(x) as.list(x)), use.names = FALSE)
        tad_parallel.end_time <- Sys.time()
        print(paste0("Parametros recobrados del chr",chr," en:"))
        (tad_parallel.time_taken <- round(tad_parallel.end_time - tad_parallel.start_time,4))
        print(paste0("Fuerza de interaccion en anclas de tads para chr",chr," calculadas"))
        return(list(mtz.chr=tads.statistics.dt))
      }
      
      list.chr <- conteos()
      
      domains.bed <- rbind(domains.bed,list.chr$mtz.chr)
      
    }else{
      print(paste0("No hay regiones en el chr",chr," que recuperar"))
    }
  } 
  
  domains.bedpe <- data.table(domains.bed[,1],
                              floor(domains.bed[,2]-(resol/2)),
                              floor(domains.bed[,2]+(resol/2)),
                              domains.bed[,1],
                              floor(domains.bed[,3]-(resol/2)),
                              floor(domains.bed[,3]+(resol/2)))
  
  colnames(domains.bedpe) <- c("chr1","start1","end1","chr2","start2","end2")
  
  domains.bedpe <- cbind(domains.bedpe,domains.bed[,-c(1:3)])
  
  iS.start <- domains.bed[,.(chr = chr, 
                             start = start-(resol/2), end = start+(resol/2),
                             id=id, Boundary="Up", id_Boundary=paste0(id,"_Up"),
                             upTAD_obs=upTAD_start_obs,
                             downTAD_obs=downTAD_start_obs,
                             interTAD_obs=interTAD_start_obs,
                             iStrength_obs=iStrength_start_obs,
                             interTAD_log=interTAD_start_log,
                             iStrength_log=iStrength_start_log,
                             interTAD_oe=interTAD_start_oe,
                             iStrength_oe=iStrength_start_oe,
                             interTAD_median=interTAD_start_median,
                             iStrength_median=iStrength_start_median)]
  iS.end <- domains.bed[,.(chr = chr, 
                           start = end-(resol/2), end = end+(resol/2),
                           id=id, Boundary="Down", id_Boundary=paste0(id,"_Down"),
                           upTAD_obs=upTAD_end_obs,
                           downTAD_obs=downTAD_end_obs,
                           interTAD_obs=interTAD_end_obs,
                           iStrength_obs=iStrength_end_obs,
                           interTAD_log=interTAD_end_log,
                           iStrength_log=iStrength_end_log,
                           interTAD_oe=interTAD_end_oe,
                           iStrength_oe=iStrength_end_oe,
                           interTAD_median=interTAD_end_median,
                           iStrength_median=iStrength_end_median)]
  
  iS.domains <- rbind(iS.start,iS.end)
  iS.domains <- iS.domains[order(chr,start)]
  iS.domains <- iS.domains[!duplicated(iS.domains[,.(chr,start)])]
  
  write.table(domains.bed,
              paste0(dir,"/iStrength/",nomen,"/",resol,"/",nomen,".iS_",window_label,".res_",resol/1000,"kb.bed"),
              row.names=FALSE,col.names=TRUE,sep="\t", quote = FALSE)
  write.table(iS.domains,
              paste0(dir,"/iStrength/",nomen,"/",resol,"/",nomen,".iS_",window_label,".res_",resol/1000,"kb.iStrength"),
              row.names=FALSE,col.names=TRUE,sep="\t", quote = FALSE)
  write.table(domains.bedpe,
              paste0(dir,"/iStrength/",nomen,"/",resol,"/",nomen,".iS_",window_label,".res_",resol/1000,"kb.bedpe"),
              row.names=FALSE,col.names=TRUE,sep="\t", quote = FALSE)
  
}