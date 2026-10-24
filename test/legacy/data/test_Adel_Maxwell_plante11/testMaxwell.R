#
#
#           Script for generating reference simulation for reference plant (Maxwell plant 11)
#            *.R files are from February10, 2013
#
#
#
#Load R files for Adel.R
#
AdelRfiles <- c("Adel.R","setAdel.R","UseAdel.R")
Adel_source_path <- file.path(
    "..", "..", "..", "src", "openalea", "adel"
)
sapply(
    file.path(Adel_source_path, AdelRfiles),
    source
)
#
#run the model with csv files
#
#read parameters from csv
manip <- c("_plante11")
chem <- c("./Maxwell_")
pars <- devTcsv(paste(chem,"axeT",manip,".csv",sep=""),
                paste(chem,"dimT",manip,".csv",sep=""),
                paste(chem,"phenT",manip,".csv",sep=""),
                paste(chem,"earT.csv",sep=""))
#Specify geoLeaf & geoAxe
geoLeaf <- genGeoLeaf()
geoAxe <- genGeoAxe()
#
attach(paste(chem,"LeafCurv_global_ntop.RData",sep=""))
xydb <- get("xy")
detach()
attach(paste(chem,"Leaf2D_global_ntop.RData",sep=""))
srdb <- get("sr")
detach()
#
#generate a list of plant to simulate from parameters
pl <- setAdel(pars$axeT, pars$dimT, pars$phenT, pars$earT, pars$ssisenT, geoLeaf, geoAxe, nplants=1, xy_db=xydb, sr_db=srdb, seed=1)
ref <- sapply(seq(0,2200,100),function(x) runAdel(x,pl)[[1]],simplify=FALSE)
write.csv(do.call('rbind',ref),'Maxwell_reference_simulation.csv',row.names=FALSE)

#run the model as a whole from plant list to AleaChn
canopy <- runAdel(-61,pl)[[1]]#
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
kinExt <- lapply(pl,function(plant) kinL(1032,plant))
kinlist <- kinExt[[1]]
d <- 1
a <- 1
kin <- data.frame(kinlist[[a]][d,,])
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
