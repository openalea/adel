library(RODBC)

#connexion fichier excel
readExcell <- function(file,onglet="Feuil1"){
res <- NULL
if (require(RODBC)) {
ch <- odbcConnectExcel(file)
res <- sqlQuery(ch,paste("select * from [",onglet,"$]",sep=""))
odbcClose(ch)}
res
}

#Precaclcul des parametres pour graphtal (die,dapp,dfeu) en fonction des rytmes d'apparition de col et de pointes
calckinfeu <- function(dev,sen,Nf) {

  smlig<-smooth.spline(dev$TT,dev$nblig)
  smvis<-smooth.spline(dev$TT,dev$nbvis)
  ddtipcol<-cbind(dev$TT,cbind(predict(smvis,dev$TT)$y, predict(smlig,dev$TT)$y))
  #ajout des points extreme
  #calcul dernier point
  dat <- ddtipcol[ddtipcol[,2]>(Nf/2) & ddtipcol[,2]<=(Nf-1),1:2]  #recup point de nbvis nf/2 à nf-1
  reg <- lsfit(dat[,2],dat[,1])  #moindre carrés y=ax+b
  tipNf <- reg$coef[1]+reg$coef[2]*Nf #TT de la derniere feuille Nf selon cette droite
  dat <- ddtipcol[ddtipcol[,3]>(Nf/2) & ddtipcol[,3]<=(Nf-1),c(1,3)] #pareil pour les cols
  reg <- lsfit(dat[,2],dat[,1])
  colNf <- reg$coef[1]+reg$coef[2]*Nf
  dat <- ddtipcol[ddtipcol[,3]<(Nf/2),c(1,3)] #recup point de nblig < nf/2
  reg <- lsfit(dat[,2],dat[,1])
  col0 <- reg$coef[1]   #donne t0 pour calcul Dil
  #
  #ajout date moyenne fin de senescence (!extrapolation après la floraison)
  smsen<-smooth.spline(sen$TT,sen$nb_sen)
  senEXTRA<- predict(smsen,seq(50,3000,100))

  res <- matrix(nrow=Nf,ncol=5)
  #Tini : date absolue d'initiation
  #die : delai initiation-extension rapide
  #dia : delai initiation - apparition tip
  #dil : duree initiation - ligulation (fin extension)
  colnames(res) <- c("Tini","Die","Dia","Dil","Dsen")
  res[,"Tini"]=approx(c(1,ddtipcol[,2],Nf),c(0,ddtipcol[,1],tipNf),c(1,(2:Nf)/2),ties="min")$y   #approx renvoie interpolations régulières
  res[,"Die"]=0.8*approx(c(1,ddtipcol[,2],Nf),c(0,ddtipcol[,1],tipNf),1:Nf,ties="min")$y-res[,"Tini"]
  res[,"Dia"]=approx(c(1,ddtipcol[,2],Nf),c(0,ddtipcol[,1],tipNf),1:Nf,ties="min")$y-res[,"Tini"]
  res[,"Dil"]=approx(c(0,ddtipcol[,3],Nf),c(col0,ddtipcol[,1],colNf),1:Nf,ties="min")$y-res[,"Tini"]
  res[,"Dsen"]= approx(senEXTRA$y,senEXTRA$x,1:Nf,ties="min")$y-res[,"Tini"]
res
}



#ecriture nbr de phytomere et check concordande entre fichiers de donnees  pour nbr de feuilles
writenumf <- function(fname,dev,morpho,anglesPrev) {
  outfile <- paste("./temp/",fname,sep="")

  nb_phy <- round(dev$nbvis[length(dev$nbvis)])
  nb_phy_morpho <- length(morpho$moy_EN)
  nb_phy_prev  <-  anglesPrev$numf[length(anglesPrev$numf)]

  if ((nb_phy == nb_phy_morpho) &  (nb_phy <= nb_phy_prev))  cat("#define NB_PHY ",nb_phy,"\n\n",file=outfile) ;   res <- "NB_PHY OK"

  if (nb_phy != nb_phy_morpho) res <- "The number of phytomers do not match between morphological and developemental measurements"

  if (nb_phy > nb_phy_prev) res = "The number of phytomers do not match between morphological and geometrical (prevot) measurements"

  res
  }



# ecriture PO
writePO <- function(fname) {
  outfile <- paste("./temp/",fname,sep="")

  cat("#define SOIL_REFLECTANCE 0.17" ,"\n\n",file=outfile,append=T)
  cat("#define LEAF_REFLECTANCE 0.09","\n\n",file=outfile,append=T)
  cat("#define LEAF_TRANSMITANCE 0.04","\n\n",file=outfile,append=T)
  cat("#define SEN_LEAF_REFLECTANCE 0.09","\n\n",file=outfile,append=T)
  cat("#define SEN_LEAF_TRANSMITANCE 0.04","\n\n",file=outfile,append=T)
  cat("#define STEM_REFLECTANCE 0.09" ,"\n\n",file=outfile,append=T)
}



# ecriture carto
writecarto <- function(fname, plan, interR, interP) {
  outfile <- paste("./temp/",fname,sep="")
  mat<-as.matrix(plan)
  dims<-dim(mat)[2]
   cat("#define NB_RANG ",dim(mat)[2],"\n\n",file=outfile,append=T)
   cat("#define NB_PLANTES ",dim(mat)[1],"\n\n",file=outfile,append=T)
   cat("#define INTER_RANG ", interR,"\n\n",file=outfile,append=T)
   cat("#define INTER_PLANTE ",interP ,"\n\n",file=outfile,append=T)
   cat("#define CARTOVALUES ",paste(t(mat),collapse=","),"\n\n",file=outfile,append=T)

}




#ecriture dim et kin et verifie type des donnes (erreur avec case vide ou NA dans les fichiers excel)
writedimkin <- function(fname, dev, sen, morpho, prec=2)  {
  outfile <- paste("./temp/",fname,sep="")
  nb_phy <- round(dev$nbvis[length(dev$nbvis)])

  if (is.numeric(morpho$moy_EN) & is.numeric(dev$nbvis)){
  txt <- "DIMKIN OK"
  #calckinfeu dans CalcPars.R
  res <- round(calckinfeu(dev,sen,nb_phy),0)
    cat("#define KINFEUVALUES ",paste(t(res),collapse=","),"\n\n",file=outfile,append=T)

  #dimensions feuilles et ecart types
  mat<-round(as.matrix(morpho[,c(1,2,4,5,7,8,10,11)]),prec-1)
  mat<-mat*0.001    #conversion en m
   cat("#define DIMVALUES ",paste(t(mat),collapse=","),"\n\n",file=outfile,append=T)
  }
  if (!is.numeric(morpho$moy_EN)) txt <- "Format du tableau de dimensions finales non numerique"

  if (!is.numeric(dev$nbvis)) txt <- "Format du tableau de developpement non numerique"

  txt
}




#ecriture prevot
writeprev<- function(fname, alpha, anglesPrev, dev, prec=2) {
  outfile <- fname#paste("./temp/",fname,sep="")
  nb_phy <- round(dev$nbvis[length(dev$nbvis)])
  nb_phy_prevMax  <-  anglesPrev$numf[length(anglesPrev$numf)]
  nb_phy_prevMin  <-  anglesPrev$numf[1]

  # valeur du facteur de forme moyen (alpha de la feuille de l'épi)
  alpha_moy<-alpha$a*2*nb_phy/3+alpha$b
  cat("#define ALPHA_LARG ",alpha_moy,"\n\n",file=outfile,append=T)

  mat <- as.matrix(anglesPrev[,c(-1:-3,-10)])
  mat <- ifelse(is.na(mat),-999,mat)
  res <- round(mat[,-1],prec)
  Nf <- nb_phy#max(mat[,"numf"])

  nblg <- aggregate(mat[,1],by=list(mat[,1]),length)$x #nbr de feuilles par rang (pour définir dim du tableau)
  #si manque la/les premiere feuilles -> duplique premier dispo pour combler les manque
  if (nb_phy_prevMin>1) for (i in 1:(nb_phy_prevMin-1)) {nblg <- c(nblg[1],nblg)}
  #si plus de feuilles que nbr median, enleve les derniere
  if (length(nblg)>nb_phy) nblg <- nblg[1:nb_phy]

  maxlg <- max(nblg) #nbval pour avoir une matrice complete (a trou !)

  cat("#define NBLGMAXGEOM ",maxlg,"\n\n",file=outfile,append=T)
  cat("#define NBLGGEOM ",paste(nblg,collapse=","),"\n\n",file=outfile,append=T)  #paste,collapse remplace séparateur ' ' pas ','
  cat("#define GEOMVALUES ",file=outfile,append=T)


  #si feuilles manquantes au debut du profil recupere premier dispo
  if (nb_phy_prevMin>1) for (i in 1:(nb_phy_prevMin-1)) {
    cat("{",file=outfile,append=T)
    dat <- res[mat[,"numf"]==nb_phy_prevMin,] #récupération de feuille 2 ds la matrice res
    cat(paste(t(rbind(dat,matrix(-999,ncol=ncol(dat),nrow=maxlg-nrow(dat)))),collapse=","),file=outfile,append=T)
    cat("}\\\n",file=outfile,append=T)       #\ de fin supplémentaire pour retour à la ligne ds les macros (spécifique graphtal)
  }

  #reste des feuilles
  for (i in nb_phy_prevMin:nb_phy) {
    if (i==1) cat("{",file=outfile,append=T)

    if (i>1) cat(",{",file=outfile,append=T)
    dat <- res[mat[,"numf"]==i,]
    cat(paste(t(rbind(dat,matrix(-999,ncol=ncol(dat),nrow=maxlg-nrow(dat)))),collapse=","),file=outfile,append=T)
    cat("}\\\n",file=outfile,append=T)
  }
  cat("\n",file=outfile,append=T)

}



#fname <- "_05_geomF2_S1B1.h"
#writenumf (fname,dev,morpho,anglesPrev)
#writePO (fname)
#writecarto(fname, plan, 0.8, 0.125)
#writedimkin(fname, dev, sen, morpho)
#writeprev(fname, alpha, anglesPrev, dev, prec=2)

