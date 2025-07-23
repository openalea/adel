#
#
#           Script for using/debugging Adel.R from within R
#
#
# Will create .RData in the directory. Don't commit it !
#
#
#
#Load R files for Adel.R
#
AdelRfiles <- c("Adel.R","setAdel.R","UseAdel.R","genString.R")
sapply(AdelRfiles,source)
#
#run the model with csv files
#
#read parameters from csv
#devT <- devTcsv('axeT.csv','dimT.csv','phenT.csv')
manip <- c("Ca0N")
chem <- c("./data/")
devT <- devTcsv(paste(chem,"axeT",manip,".csv",sep=""),
                paste(chem,"dimT",manip,".csv",sep=""),
                paste(chem,"phenT",manip,".csv",sep=""),
                paste(chem,"earT",manip,".csv",sep=""),
                paste(chem,"ssi2sen.csv",sep=""))
#
#avatar
#devT <- devTcsv(paste(chem,"axeTSoissonsNormal_10plants.csv",sep=""),
#                paste(chem,"dimTSoissons.csv",sep=""),
#                paste(chem,"phenTNormalPhyllochron.csv",sep=""),
#                paste(chem,"earT",manip,".csv",sep=""),
#                paste(chem,"ssi2sen.csv",sep=""))
#Specify geoLeaf & geoAxe
geoLeaf <- genGeoLeaf()
geoAxe <- genGeoAxe()
#
attach(paste(chem,"So99.RData",sep=""))
xydb <- get("So99")
detach()
attach(paste(chem,"SRSo.RData",sep=""))
srdb <- get("SRSo")
detach()
#
#generate a list of plant to simulate from parameters
pl <- setAdel(devT$axeT,devT$dimT,devT$phenT,devT$earT,devT$ssisenT,geoLeaf,geoAxe,nplants=1,xy_db=xydb,sr_db=srdb)
#
# load('adel_pars.RData')
#pl <- plants
#run the model as a whole from plant list to AleaChn
canopy <- runAdel(1200,pl)[[1]]
#
Scanopy <- canL2canS(canopy,srdb)
#
chn <- genString(canopy)
#
#
#Run step by step the model
#
#
#Predict kinetic of extension
#
kinExt <- lapply(pl,function(plant) kinL(2000,plant))
#kinExt[[1]][[1]][1,,]
#
# idem for visible parts
kinExtv <- lapply(kinExt,kinLvis)
#
#Dataframe representing the canopy
#
canopy <- getdesc(kinExtv,pl)
#
#string
#
chn <- genString(canopy)
#
#
#Test setAdelArv
#
load("./data/Arvalis.RData")
source("./data/Rfun.R")
source("ArvalisToAdel.R")
#
# simulates how adel param is relaed to HS in case collar app = endleaf
#
phyl <- 110
endLeaf <- 1.6
nF <- 11
#
time <- seq(0,1500) # time since emergence leaf 1
#
tips <- 1 + time / phyl
tips[1] = 0
#
colars = 1 +time / phyl - endLeaf
colars[time <= endLeaf *phyl] = time[time <= endLeaf *phyl] / (endLeaf * phyl)
leaves <- matrix(0,ncol=nF,nrow=length(time))
for (i in seq(nF)) {
  xem <- phyl * (i - 1)
  xlig <- xem + endLeaf * phyl
  l <- leaves[,i]
  l[time >= xem] <- (time[time >= xem] - xem) / phyl / endLeaf
  l[time >= xlig] <- 1
  leaves[,i] <- l
}

#
plot(c(-1,7), c(-1,6),xlab='phyllochronic time since plant emergence', ylab='Leaf #',pch='',cex.lab=1.5)
abline(0,1,lty=3,col=1)
#leaf start
abline(1.4,1, lty=2, col=5)
# tips
lines(time/phyl, tips, col=4, lwd=2)
#colars
lines(time/phyl, colars, col=3, lwd=2)
#tillers
abline(-2,1,col=6)
#
# Klepper's shoot events
for (event in list(list(0,1,'L1 appears',4), list(1,2,'L2 starts',4), list(2.7,1,'T1 starts',2), list(3,1,'T1 appears',4))) {
x = event[1]
y = event[2]
t = event[3]
p = event[4]
points(x,y,pch=20,col=2)
text(x,y,t, cex=0.7,col=2, pos=p)
}
#
# Klepper's root events
for (event in list(list(-1,-1,'-2AB start',4),list(0.2,-1,'-1AB start',4), list(1.6,0,'0AB start',4),list(2.8,1,'1AB start',3), list(4.3,1,'1x starts',4), list(5.3,1,'1AB Branch',4), list(6.8,1,'1x Branch',1))) {
x = event[1]
y = event[2]
t = event[3]
p = event[4]
points(x,y,pch=20,col=4)
text(x,y,t, cex=0.7,col=4, pos=p)
}
#
# base legend
legend(-1,6,c('Tips','Collars (Fournier et al.)','Leaf starts (Fournier et al.)','Tillers (Kirby et al)'),col=c(4,3,5,6),lty=c(1,1,2,1),lwd=c(2,2,1,1))
#
# adding models and update legend
# ABroot: collars n+1
lines(time/phyl, colars -1, col=3, lty=2)
# AB root branches : 2.5 phyllo after AB starts
lines(time/phyl+2.5, colars -1, col=3, lty=3)
# child tiller colars : 1.6 phyllocron after tiller emerges (x root starts on parent axis)
lines(time/phyl + 2.6, colars, col=6, lty=2)
# 1x branches
lines(time/phyl+5.1, colars, col=6, lty=3)
legend(-1,6,c('Tips','Collars (Fournier et al.)','Leaf starts (Fournier et al.)','Tillers (Kirby et al)', 'Tiller + AB start (this study)', 'x start (this study)'),col=c(4,3,5,6,3,6),lty=c(1,1,2,1,2,2),lwd=c(2,2,1,1,1,1))
#
# adding HaunStage 
#
hs <- apply(leaves,1,function(x) {ifelse(any(x>=1),max(which(x>=1)) + sum(x[(max(which(x>=1))+1):length(x)]),sum(x))})
hst <- hs
sel <- hs < 1
hst[sel] <- (time[sel]/phyl + 2.6 - 3) / 1.6
hst[time/phyl + 2.6 - 3< 0] <- 0 
#
#TODO : arrange the fact that HStilller = zero at tiller appearance (1st HS is 0.4 pyllochron shorter on tillers !)
lines(time / phyl,hs, lwd=2, col=2)
lines(time / phyl + 2.6,hst, lwd=2, col=6)
lines(time / phyl + 2.6,hs, lty=3, col=6)
tp = (time / phyl)
sel = tp > endLeaf & tp < (nF - endLeaf)
fit <- lsfit(tp[sel], hs[sel])
abline(fit$coef[1],fit$coef[2],lty=3,col=2)
text(4,-0.5,paste('delay Tip->HS:', -round(mean(hs[sel] - tp[sel] - 1),2)))
text(4,-1,paste('first HS of child tiller last 1.2 phylllo instead of 1.6', -round(mean(hs[sel] - tp[sel] - 1),2)))
legend(-1,6,c('HS','Tips','Collars','leaf start','child tiller tip','child tiller HS'),col=c(2,4,3,5,6,6),lty=c(1,1,1,2,1,1),lwd=c(2,2,2,1,1,2))
#
# Expectation
#
# Crown roots = f(HS)
#
ab = colars - 1
ab[colars < 2] = 0
plot(hs,2*ab,xlab='Haun stage', ylab='Root #',pch='',cex.lab=1.5) 
lines(hs,2*ab, col=2)
abline(-1.57*1.95,1.95, col=3)
legend(1,20,c('predicted', 'Klepper'),col=c(2,3), lty=1)