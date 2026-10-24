#
# User front End and utilities for using adel
#
#
#
#Run adel for several dates and a list of plants. Returns a list of string
#
runAdel <- function(dates,plants,pars = list('senescence_leaf_shrink' = 0.5, 'leafDuration' = 2, 'fracLeaf' = 0.2, 'stemDuration' = 2 / 1.2, 'dHS_col'=0.2, 'dHS_en'=0.,'epsillon' = 1e-6, 'HSstart_inclination_tiller' = 1, 'rate_inclination_tiller' = 30,'drop_empty'=TRUE)) {
  out <- vector("list",length(dates))
  #deals with python-flatten lists
  if ("axeT" %in% names(plants)) {
    plants <- list(plants)
    names(plants) <- plants[[1]]$refp
  }
  for (i in seq(out)) {
    kinlist <- lapply(plants,function(plant) kinLvis(kinL(dates[i],plant,pars),pars))
    desc <- getdesc(kinlist,plants,pars)
    #chn <- genString(desc,pars)
    if (!is.null(desc))
      out[[i]] <- cbind(TT=dates[i],desc)
  }
  out
}
#
#set Adel from parameters 
#
setAdeluser <- function(devT,geoLeaf,geoAxe,nplants,sample='random',seed=NULL,xy_db=NULL,sr_db=NULL, ssipars=NULL) {
  setAdel(devT$axeT,devT$dimT,devT$phenT,devT$earT,devT$ssisenT,geoLeaf,geoAxe,nplants,sample,seed,xy_db,sr_db)
}
#
#build devT from csv parameter files
#
readCsv <- function(file,type=1) {
  reader <- get(ifelse(type==1,"read.csv","read.csv2"))
  #type 1 : "." for decimal, "," for separator
  #type 2:  "," for decimal, ";" for separator
  data <- reader(file)
  #filter empty columns
  Filter(function(x)!all(is.na(x)), data)
}
#
devTcsv <- function(axeTfile,dimTfile,phenTfile,earTfile=NULL,ssisenTfile=NULL,type=1) {

  #conversion nouvelle nomencalture
  #axeT
  axeT <- readCsv(axeTfile)
  # nouvelle convention (kirby) id_axe
  if ("axe" %in% colnames(axeT))
    axeT$axe <- ifelse(axeT$axe==0,"MS",paste("T",axeT$axe,sep=""))
  #
  #conversions nouveaux noms
    
  conv <- c("plant","axe","nf","end","disp","dimIndex","phenIndex","earIndex","emf1","ligf1","senf1","dispf1")
  names(conv) <- c("id_plt","id_axis","N_phytomer","TT_stop_axis","TT_del_axis","id_dim","id_phen","id_ear","TT_em_phytomer1","TT_col_phytomer1","TT_sen_phytomer1","TT_del_phytomer1")
  if (all(conv %in% colnames(axeT)))
    colnames(axeT)[colnames(axeT) %in% conv] <- names(conv)[na.omit(match(colnames(axeT),conv))]
  else if (!all(names(conv) %in% colnames(axeT)))
    stop(paste("axeT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(axeT)],collapse=" ")))
  # force numeric conversion
  numcols <- names(conv)[-grep('id_axis',names(conv))]
  for (w in numcols)
    axeT[,w] <- as.numeric(as.character(axeT[,w]))
  #dimT
  dimT <- readCsv(dimTfile, type)
  conv <- c("index","nrel","Ll","Lw","Gl","Gd","El","Ed")
  names(conv) <- c("id_dim","index_rel_phytomer","L_blade","W_blade","L_sheath","W_sheath","L_internode","W_internode")
  if (all(conv %in% colnames(dimT)))
    colnames(dimT)[colnames(dimT) %in% conv] <- names(conv)[na.omit(match(colnames(dimT),conv))]
  else if (!all(names(conv[-grep('nrel',conv)]) %in% colnames(dimT)))
    stop(paste("dimT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(dimT)],collapse=" ")))
  # phenT
  phenT = readCsv(phenTfile, type)
  conv <- c("index","nrel","tip","col","ssi","disp")
  names(conv) <- c("id_phen","index_rel_phytomer","dTT_em_phytomer","dTT_col_phytomer","dTT_sen_phytomer","dTT_del_phytomer")
  if (all(conv %in% colnames(phenT)))
    colnames(phenT)[colnames(phenT) %in% conv] <- names(conv)[na.omit(match(colnames(phenT),conv))]
  else if (!all(names(conv[-grep('nrel',conv)]) %in% colnames(phenT)))
    stop(paste("phenT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(phenT)],collapse=" ")))
  phenT <- phenT[!is.na(phenT$id_phen),]
  # earT
  if (!is.null(earTfile)) {
    earT <- readCsv(earTfile, type)
    conv <- c("index","em_ear","em_ped","end_gf","l_ped","d_ped","l_ear","Sp_ear","l_ear_awn")
    names(conv) <- c("id_ear","dTT_em_ear","dTT_em_peduncle","TT_z92","L_peduncle","W_peduncle","L_ear","A_ear","L_spike")
    if (all(conv %in% colnames(earT)))
    colnames(earT)[colnames(earT) %in% conv] <- names(conv)[na.omit(match(colnames(earT),conv))]
  else if (!all(names(conv) %in% colnames(earT)))
    stop(paste("earT : missing data: ",paste(names(conv)[!names(conv) %in% colnames(earT)],collapse=" ")))
  }
  else
    earT <- NULL
    
  if (!is.null(ssisenTfile))
    ssisenT <- readCsv(ssisenTfile, type)
  else
    ssisenT <- NULL
  list(axeT = axeT,
       dimT = dimT,
       phenT = phenT,
       earT = earT,
       ssisenT = ssisenT)
}
#
#
#geoAxe from parameter
#
genGeoAxe <- function(azTM = 75,dazT = 5,incBmM = 2,dincBm = 2,incT = 60,dincT = 5,depMax = 7,depMin = 1.5,dazTb=60) {
  list(
       azT = function(a) {
         ifelse(a == 'MS',
                runif(1) * 360,#plant azimuth
                azTM + (runif(1) - .5) * dazT)
       },
       azTb = function(a) {
         ifelse(a == 'MS',
                0,
                (runif(1) - .5) * dazTb)
       },
       incT = function(a) {
         ifelse(a == 'MS',
                incBmM + (runif(1) - .5) * dincBm,
                incT + (runif(1) - .5) * dincT)
       },
       dredT = function(a) {
         #1.5 is an offset to avoid tiller superposed to mainstem
         ifelse(a == 'MS',
                0,
                depMin + runif(1) * (depMax-depMin))
       }
       )
}
#
#geoLeaf from parameter (prevoir aussi une boite freeGeomAxe et freegeoLeaf)
#
genGeoLeaf <- function(ntoplim = 4,dazTop = 60,dazBase = 30,topIndex=TRUE) {
   list(
        Azim = function(a,n,nf) {
          ntop = nf - n
          ifelse(ntop <= ntoplim,
                 180 + dazTop * (runif(1) - .5),
                 180 + dazBase * (runif(1) - .5))
               },
        Lindex = function(a,n,nf) {
          ifelse(topIndex,
                 nf - n + 1,
                 n)}
        )
 }
#
# Compute surfaces from lengths in canopy table
#
leafSurface <- function(shape_db, shape_index, scL,scW,from=0,to=1) {
  shape <- shape_db[[shape_index]]
  if (is.na(scL) | is.na(scW) | is.null(shape))
    res <- NA
  else if (to <= from | scL == 0 | scW == 0)
    res <- 0
  else {
    shape[,1] <- shape[,1] / max(shape[,1]) * scL
    shape[,2] <- shape[,2] / max(shape[,2]) * scW
    ws <- approxfun(shape[,1],shape[,2],rule=2)
    res <- integrate(ws,from * scL,to *scL)$value
  }
  res
}
#
canL2canS <- function(canL,sr_db,leaf_shrink=NULL) {
  canS <- canL
  canS$ntop <- canL$nff - canL$numphy + 1
  if (all(canS$LcIndex <= 1)) {#leaf shape has not been set by setAdel
    print("canL2canS : can't compute Blade surfaces : SR data not connected to setAdel!!!")
  }
  else {
    rg <- ifelse(is.na(canS$LcType) | canS$LcType <= 0, 1,canS$LcType)#handle phytomers above nff (peduncle, ear, awn) for which Lindex = zero
    shapes <- sr_db
    Lref <- canS$L_shape
    Lwref <- canS$Lw_shape
    Lwsen <- canS$Lw_shape * canS$LsenShrink
    # add some info on visibility
    canS$Lvsen <- pmin(canS$Lv,canS$Lsen)
    canS$Lvgreen <- canS$Lv - canS$Lvsen
    canS$Gvsen <- pmin(canS$Gv,canS$Gsen)
    canS$Gvgreen <- canS$Gv - canS$Gvsen
    canS$Evsen <- pmin(canS$Ev,canS$Esen)
    canS$Evgreen <- canS$Ev - canS$Evsen
    #
    #add info on distances base of axes -> col
    g <- list(canS$plant,canS$axe_id)
    canS <- unsplit(lapply(split(canS,g), function(dat) {dat$d_basecol <- cumsum(dat$El) + dat$Gl;dat}),g)
    #
    canS$S_shape <- sapply(seq(nrow(canS)),function(x) leafSurface(shapes, rg[x], Lref[x], Lwref[x]))
    #
    base <- 1 - (canS$Lvsen + canS$Lvgreen) / Lref
    top <- 1 - canS$Lvsen / Lref
    canS$Slvgreen <- sapply(seq(nrow(canS)),function(x) leafSurface(shapes, rg[x], Lref[x], Lwref[x], base[x], top[x]))
    #
    base <- 1 - canS$Lvsen / Lref
    canS$Slvsen <- sapply(seq(nrow(canS)),function(x) leafSurface(shapes, rg[x], Lref[x], Lwsen[x], base[x], 1))
    #
    canS$Slv <- canS$Slvgreen + canS$Slvsen
    #
    top <- 1 - canS$Lsen / Lref
    Shgreen <- sapply(seq(nrow(canS)),function(x) leafSurface(shapes, rg[x], Lref[x], Lwref[x], 0, top[x]))
    #
    base <- 1 - canS$Lsen / Lref
    canS$Slsen <- sapply(seq(nrow(canS)),function(x) leafSurface(shapes, rg[x], Lref[x], Lwsen[x], base[x], 1))
    #
    canS$SLl <- Shgreen + canS$Slsen    
  }
  #names(res)[match(c("Ll","Lsen","Lv"),names(res))] <- c("SLl","SLsen","SLv")
  #
  for (w in c("Gl","Gv","Gsen","Gvsen","Gvgreen"))
    canS[[paste("S",w,sep="")]] <- canS[[w]] * pi * canS$Gd
  for (w in c("El","Ev","Esen","Evsen","Evgreen"))
    canS[[paste("S",w,sep="")]] <- canS[[w]] * pi * canS$Ed
  canS
}
             
                         
    
    
