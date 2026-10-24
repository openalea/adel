#
#
#          R routine for converting a dataframe representation of canopy to String for CanMTG interpreter
#
#
# conventions tissue_ttype :
ttype <- list(green_lamina=1,sen_lamina=2,green_sheath=3,sen_sheath=4,green_internode=5, sen_internode=6,green_peduncle = 7, sen_peduncule=8, green_ear=9, sen_ear=10, green_awn = 11, sen_awn = 12)
# 1 = Green lamina
# 2 = Senescent Lamina
# 3 = Green sheath
# 4 = Senescent sheath1,3,
# 5 = Green internode
# 6 = senescent Internode
# 7 = green peduncule
# 8 = senescent peduncule
# 9 = green ear
# 10 = senescent ear
# 11 = green awn
# 12 = senescent awn
StemElement <- function(po, length, dbase, dtop,prec=3) {
  paste("StemElement(",paste(c(po,round(c(length,dbase,dtop),prec)),collapse=","),")",sep="")
}
#
LeafElement <- function(po, Lf, l, wM, srb, srtop,index,seed,incB,prec=3) {
  paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed,incB),prec)),collapse=","),")",sep="")
  #paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed),prec)),collapse=","),")",sep="") # incB removed untill debugging completed
}
indexedLeafElement <- function(po, Lf, l, wM, srb, srtop,index,seed,incB,tindex,prec=3) {
  paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed,incB),prec),tindex),collapse=","),")",sep="")
  #paste("LeafElement(",paste(c(po,round(c(Lf,l,wM,srb,srtop),prec),index,round(c(seed),prec)),collapse=","),")",sep="") # incB removed untill debugging completed
}
#
pitch <- function(inc) paste("+(",round(inc,2),")",sep="")
roll <- function(ang) paste("/(",round(ang,2),")",sep="")
up <- function(ang) paste("^(",round(ang,2),")",sep="")
#
Internode <- function(lv,lvsen,diam,po,pos,epsillon) {
  chn <- ""
  lvgreen = lv - lvsen
#  if (lvgreen > epsillon) # no condition, as internode should always be there to allow bearing tillers
    chn <- paste(chn,StemElement(ttype$green_internode,lvgreen,diam,diam),sep="")
  #if (lvsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_internode,lvsen,diam,diam))
  chn
}
#
Peduncle <- function(lv,lvsen,diam,po,pos,epsillon) {
  chn <- ""
  lvgreen = lv - lvsen
  if (lvgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_peduncle,lvgreen,diam,diam),sep="")
  if (lvsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_peduncle,lvsen,diam,diam))
  chn
}
#
Ear <- function(lv,lvsen,diam,po,pos,epsillon) {
  chn <- ""
  lvgreen = lv - lvsen
  if (lvgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_ear,lvgreen,diam,diam),sep="")
  if (lvsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_ear,lvsen,diam,diam))
  chn
}
#
Awn <- function(lv,lvsen,diam,po,pos,epsillon) {
  chn <- ""
  lvgreen = lv - lvsen
  if (lvgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_awn,lvgreen,diam,diam),sep="")
  if (lvsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_awn,lvsen,diam,diam))
  chn
}
#
Sheath <- function(lv,lvsen,diam,po,pos,epsillon) {
  chn <- ""
  lvgreen = lv - lvsen
  if (lvgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_sheath,lvgreen,diam,diam),sep="")
  if (lvsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_sheath,lvsen,diam,diam))
  chn
}
#
RolledBlade <- function(lv,lvsen,diam,po,pos,epsillon) {
  chn <- ""
  lvgreen = lv - lvsen
  if (lvgreen > epsillon)
    chn <- paste(chn,StemElement(ttype$green_lamina,lvgreen,diam,diam),sep="")
  if (lvsen > epsillon)
    chn <- paste(chn,StemElement(ttype$sen_lamina,lvsen,diam,diam))
  chn
}
#
Blade <- function(lv, lvsen, Lf, wM, lsenshrink,lctype, lcindex, incB, po, pos, epsillon,randindex) {
  chn <- ""
  lvgreen = lv - lvsen
  if (randindex) {#lcindex has not been set by adel)
    if (lvgreen > epsillon)
      chn <- paste(chn,
                   LeafElement(ttype$green_lamina, Lf, lv, wM, 0, lvgreen / lv,lctype, lcindex, incB))
    if (lvsen > epsillon)
      chn <- paste(chn,
                   LeafElement(ttype$sen_lamina, Lf, lv, wM * lsenshrink, lvgreen / lv, 1, lctype, lcindex, incB))
  } else {
    if (lvgreen > epsillon)
      chn <- paste(chn,
                   indexedLeafElement(ttype$green_lamina, Lf, lv, wM, 0, lvgreen / lv,lctype, lcindex, incB,lcindex))
    if (lvsen > epsillon)
      chn <- paste(chn,
                   indexedLeafElement(ttype$sen_lamina, Lf, lv, wM * lsenshrink, lvgreen / lv, 1, lctype, lcindex, incB,lcindex))
  }
 
  chn
}
#
Metamer <- function(dat,epsillon,azcum,axil = NULL) {
  
  azm <- (azcum + dat$Laz) %% 360
  
  chn <- "newMetamer"
  if (abs(dat$Einc) > 0)
    chn <- paste(chn,
                 up(dat$Einc))

  #if (dat$Ev > epsillon) #not tested to ensure persistence of stems even with 0 length (for tiller attachment)
  chn <- paste(chn,
                 Internode(dat$Ev,min(dat$Esen,dat$Ev),dat$Ed,dat$Epo,dat$Epos,epsillon))

  if (!is.null(axil))
    for (i in seq(along=unique(axil$axe_id))) {
      ax <- axil[axil$axe_id == unique(axil$axe_id)[i],]
      azaxil <- azm + ax$az_insertion[1]
      #if (i > 1)
      #  azaxil = azaxil + 170
      chn <- paste(chn,
                 "[",
                 roll(azaxil),
                 "newAxe")
      nf <- ax$nff[1]
      for (n in ax$numphy[1:min(nrow(ax),nf)]) {
        #print(n)
        chn <- paste(chn,Metamer(ax[ax$numphy == n,],epsillon,azaxil))
        azaxil = azaxil + ax$Laz[ax$numphy == n]
      }
      #add Peduncle
      if (nrow(ax) > nf)                      
        if (ax$Ev[nf + 1] > epsillon)
          chn <- paste(chn,
                 Peduncle(ax$Ev[nf + 1],min(ax$Esen[nf + 1],ax$Ev[nf + 1]),ax$Ed[nf + 1],ax$Epo[nf + 1],ax$Epos[nf + 1],epsillon))
      #add ear
      if (nrow(ax) > (nf + 1))
        if (ax$Ev[nf + 2] > epsillon)
          chn <- paste(chn,
                 Ear(ax$Ev[nf + 2],min(ax$Esen[nf + 2],ax$Ev[nf + 2]),ax$Ed[nf + 2],ax$Epo[nf + 2],ax$Epos[nf + 2],epsillon))
      #add awn
      if (nrow(ax) > (nf + 2))
        if (ax$Ev[nf + 3] > epsillon)
          chn <- paste(chn,
                 Awn(ax$Ev[nf + 3],min(ax$Esen[nf + 3],ax$Ev[nf + 3]),ax$Ed[nf + 3],ax$Epo[nf + 3],ax$Epos[nf + 3],epsillon))
          
      chn <- paste(chn,
                 "]")
    } 


  if (abs(dat$Ginc) > 0)
    chn <- paste(chn,
                 up(dat$Ginc))

  if (dat$Gv > epsillon)
    chn <- paste(chn,
                 Sheath(dat$Gv,min(dat$Gsen,dat$Gv),dat$Gd,dat$Gpo,dat$Gpos,epsillon))

  if (dat$Lr > epsillon)
    chn <- paste(chn,
                 RolledBlade(dat$Lr,min(dat$Lsen,dat$Lr),dat$Gd,dat$Lpo,dat$Lpos,epsillon))
  Lvflat = dat$Lv - dat$Lr
  if (Lvflat > epsillon) 
    chn <- paste(chn,
                 "[",
                 roll(azm),
                 Blade(Lvflat, min(dat$Lsen,Lvflat), dat$Ll,dat$Lw,dat$LsenShrink,dat$LcType,dat$LcIndex,dat$Linc,dat$Lpo,dat$Lpos,epsillon,all(dat$LcIndex <= 1) & !all(dat$LcIndex == 1)),
                 "]")
  
  chn
}
#
#axe_code = 0 pour MS, num bearing metamer pour Tiller
#
axe_code <- function(axeid) {
  idcode <- strsplit(axeid,split=".",fixed=TRUE)[[1]][1]
  if (idcode=="MS")
    code <- 0
  else
    code <- max(1,as.numeric(strsplit(idcode,split="T")[[1]][2])) # T0 attached to phyto1
  code
}

#
genString <- function(can,pars=list("epsillon" = 1e-6)) {
  epsillon = pars$epsillon
  chn <- ""
  if (!is.null(can)) {#no plant is emerged
    for (p in unique(can$plant)) {
      pl <- can[can$plant == p,]
      pl$axe <- sapply(as.character(pl$axe_id),axe_code)
      axepos <- unique(pl$axe[pl$axe_id != "MS"])
      axe <- pl[pl$axe_id == 'MS',]#main stem
      if (nrow(axe) > 0) {
        azcum <- 0
        chn <- paste(chn,
                 "[ newPlant newAxe")
                                        # stringify main stem and delegates to Metamer axilary branches, if any
        nf <- axe$nff[1]
        for (m in axe$numphy[1:min(nrow(axe),nf)]) {
          #print(m)
          if (m %in% axepos)
            chn <- paste(chn,Metamer(axe[axe$numphy == m,],epsillon,azcum,pl[pl$axe == m,]))
          else
            chn <- paste(chn,Metamer(axe[axe$numphy == m,],epsillon,azcum))
          azcum <- azcum + axe$Laz[axe$numphy == m]
        }
       
        #add Peduncle
        if (nrow(axe) > nf) 
          if (axe$Ev[nf + 1] > epsillon)
            chn <- paste(chn,
                 Peduncle(axe$Ev[nf + 1],min(axe$Esen[nf + 1],axe$Ev[nf + 1]),axe$Ed[nf + 1],axe$Epo[nf + 1],axe$Epos[nf + 1],epsillon))
        #add ear
        if (nrow(axe) > (nf + 1)) 
          if (axe$Ev[nf + 2] > epsillon)
            chn <- paste(chn,
                 Ear(axe$Ev[nf + 2],min(axe$Esen[nf + 2],axe$Ev[nf + 2]),axe$Ed[nf + 2],axe$Epo[nf + 2],axe$Epos[nf + 2],epsillon))
        #add awn
        if (nrow(axe) > (nf + 2))
          if (axe$Ev[nf + 3] > epsillon)
            chn <- paste(chn,
                 Awn(axe$Ev[nf + 3],min(axe$Esen[nf + 3],axe$Ev[nf + 3]),axe$Ed[nf + 3],axe$Epo[nf + 3],axe$Epos[nf + 3],epsillon))
        
        chn <- paste(chn,"]")
      }
    }
  }
  chn
}

