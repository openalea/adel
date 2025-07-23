#
#                Codage adel simexp
#
#
#A. Utilitaires de generation de tables de parametres ADEL a partir des donnees arvalis
#
#A.1 Generateurs pour le tableau axeTable
#
# AxeTable =  1 ligne par axe avec
# * axe :position topologique de l'axe dans la plante (0 = bm,1 = portee par phyto 1...) 
# * emf1: date d'apparition axe (emergence pointe premiere feuille)
# * end: date d'arret de croisssance de l'axe (si disparition, NA sinon)
# * disp : date de disparition de l'axe (NA = ne disparait pas)
# * azT: azimuth plan de l'axe / plan de l'axe porteur (deg)
# * incT inclinaison de l'axe (a la base) / axe porteur (deg)
# * dredT distance a floraison entre epi de l'axe et epi axe porteur
#
#
# A.1.1 Fonction renvoyant les parametres geometriques d'axe apartir d'un vecteur de position topo des axes
#
geomAxeArv <- function(axes,RFun) {
  geom = sapply(axes,RFun$geomT,simplify=F)
  azT <- NULL;incT <- NULL;dred <- NULL
  for (i in seq(axes)) {
    azT <- c(azT,geom[[i]]$azT)
    incT <- c(incT,geom[[i]]$incT)
    dred <- c(dred,geom[[i]]$dred)
  }
    
  out <- data.frame(axe = axes, azT = azT,incT = incT, dredT = dred)
  out
}
#
# A.1.2 constructeur d'un generateur de parametre pheno d'axe a partir d'abondance moyenne des axes (Axedb) et du modele de tallage. 
#
buildgenAxe <- function(axedata, phyl = 110, delreg = 2, dureg = 3,deldisp = 2,posv = 0:7, decv = c(0,3.05,3.6,4.6,6,7,8),nf = round(11 - c(0,3.05,3.6,4.6,6,7,8))) {
  emis <- axedata[1,] > 0
  naxedeb <- axedata[1,emis]
  naxefin <- axedata[2,emis]
  #fractions de talles dans naxedeb = proba d'etre sur la plante pour la talle fractionaire
  pdeb <- naxedeb - floor(naxedeb)
  #proba d'etre presentes a la fin pour les talles presentes au debut = naxefin/naxedeb
  pfin <- naxefin / naxedeb
  
  #generateur
  genAxe <- function(sTem = 0) {
    #passage en nombre entier d'axe emis selon pdeb
    naxedeb <- floor(naxedeb) + ifelse(pdeb >= runif(pdeb),1,0)
    #calcul nbre de talles presente a la fin selon pfin
    naxefin <- naxedeb * pfin
    frac <- naxefin - floor(naxefin)
    naxefin <- floor(naxefin) + ifelse(frac >= runif(frac),1,0)

    #pheno de l'axe selon modele
    nv <- max(which(naxedeb > 0))
    debreg <- sTem + (decv[nv] + delreg) * phyl
    endreg <- sTem + (decv[nv] + delreg + dureg) * phyl
    vagues <- seq(nv)
    axes <- rep(vagues,naxedeb[vagues])
    
    disp <- NULL
    for (v in vagues) {
      d <- runif(naxedeb[v],debreg,endreg)
      d[seq(length=naxefin[v])] <- NA
      disp <- c(disp,d)
    }
    data.frame(axe = posv[axes],nf = nf[axes],emf1 = sTem + decv[axes] * phyl,end = disp - deldisp * phyl,disp = disp)
  }

  return(genAxe)
}
#
#
# A.2 construction d'un generateur de parametre phytomer (dimension + azimuth + shape leaf) a partir de dimensions cal�es
#
buildgenDim <- function(dimF, Rfun,indexes = 1:6) {
  dim <- dimF[,c("L","l","Lg","d","Le","d","Phi0","StemRate","StemRate"),]
  dimnames(dim)[[2]] <- c("Ll","Lw","Gl","Gd","El","Ed","Azim","Lindex","Lseed")
  dimnames(dim)[[3]] <- (0:10)[seq(dim(dimF)[3])]
  #
  dim[,c("Gd","Ed"),] <- dim[,c("Gd","Ed"),] / 10

  genDim <- function() {
    out <- dim
    for (a in seq(dim(dim)[3])) {
      nf <- length(na.omit(dim[,1,a]))
      out[seq(nf),"Azim",a] <- sapply(seq(nf),function(n) Rfun$azim(a-1,n,nf-n))
      out[seq(nf),"Lindex",a] <- sapply(seq(nf),function(n) Rfun$xysr(a-1,n,nf-n))
      out[seq(nf),"Lseed",a] <- runif(nf)
    }
    out
  }
  genDim
}
#
#A.3 construction d'un generateur de fonctions pheno (tip,col) a  partir du phyllochrone et du tableau des axes
#
buildgenPhen <- function(phyl=110,nfveg=4.5,nfflo=2.7,nfend=0,dphylflo=5.2,dttendsen=550) {

  gentipcol <- function(axeT) {
    nv <- nrow(axeT)
    
    tip <- function(x) {
      res <- matrix(NA,nrow=length(x),ncol=nv)
      colnames(res) <- axeT$axe
      for (a in seq(nv)) {
        res[,a] <- 1 + (x - axeT$emf1[a]) / phyl
        #si end non null, on stope le dvpt a end
        if (!is.na(axeT$end[a]))
          res[x > axeT$end[a],a] <-  1 + (axeT$end[a] - axeT$emf1[a]) / phyl
      }
      res
    }
    
    col <- function(x) tip(x) - 1.6

    ssi <- function(x) {
      res <- tip(x)
      for (a in seq(ncol(res)))
        res[,a] <- approx(c(nfveg*phyl,axeT$nf[a]*phyl,(axeT$nf[a]+dphylflo)*phyl,(axeT$nf[a]+dphylflo)*phyl+dttendsen),c(0,axeT$nf[a]-nfveg,axeT$nf[a]-nfflo,axeT$nf[a]),xout=x,rule=2)$y
      res
    }

    disp <- function(x) ssi(x) - 2

    list(tip=tip,col=col,ssi=ssi,disp=disp)
  }
  
  genPhen <- function(axeT) {
    tipcol <- gentipcol(axeT)
    out <- vector("list",nrow(axeT))
    names(out) <- axeT$axe
    x <- -1000:2000
    tip <- tipcol$tip(x)
    col <- tipcol$col(x)
    ssi <- tipcol$ssi(x)
    disp <- tipcol$disp(x)
    for (a in seq(out))
      out[[a]] <- data.frame(n = 0:axeT$nf[a],tip = approx(tip[,a],x,0:axeT$nf[a],rule = 2)$y,col = approx(col[,a],x,0:axeT$nf[a],rule=2)$y,ssi = approx(ssi[,a],x,0:axeT$nf[a],rule=2)$y,disp = approx(disp[,a],x,0:axeT$nf[a],rule=2)$y)
    out
  }
  genPhen
}
#
Rfuntest = list(azim=function(a,n,ntop) {0 * runif(1)})
#
# A.4 : generateur de parametres adel a partir des donn�es arvalis
#
setAdelArv <- function(calage,nplants=1,sdlevee=0,RFun) {
  phyl <- calage$phen$PHYL
  nf <- apply(calage$dim[,1,],2,function(x) length(na.omit(x)))
  genAxe <- buildgenAxe(calage$axe,phyl,nf=nf)
  genDim <- buildgenDim(calage$dim,RFun)
  genPhen <- buildgenPhen(phyl)
  out <- vector("list",nplants)
  sTem = rnorm(nplants,sd=sdlevee)
  for (p in seq(out)) {
    axeT <- genAxe(sTem=sTem[p])
    geom <- geomAxeArv(axeT$axe,RFun)
    phenT <- genPhen(axeT)
    pedT <- phenT
    for (a in seq(pedT))
      pedT[[a]] <- data.frame(startPed=max(phenT[[a]]$col), endPed = max(phenT[[a]]$col) + 100, senPed = max(phenT[[a]]$col) + 1000)
    
    out[[p]] <- list(axeT = merge(axeT,geom),phytoT=genDim(),pheno=phenT,pedT=pedT,ssisenT=data.frame(ndel=1:4,rate1=0.07,dssit1=c(0,1.2,2.5,3),dssit2=c(1.2,2.5,3.7,4)))
  }
  out
}

