#
#
#                     Description of Adel's inputs
#
#
# Adel has three kinds of inputs available for the user:
#
# * inputs charactering the development of plants in the canopy (topology, developmental rate and size of plants)
# * inputs characterising the geometry of axes and of leaves
# * inputs charactering the simulation (time step..) and the plot configuration (number of plant, position)

# Beside these user inputs, Adel also has 'fixed' parameters and methodologies (set and 'readable' in Adel.R, but not documented) describing coordination of leaves, dynamics of geometry, computation of visibility and progression of senescence as a function of ssi. See Robert et al (Functional Plant Biology 2008) for description.

# A. Inputs describing development
# -------------------------------

#Development is parameterised in a Rlist containing 3 tables (axeT,dimT and phenT) . Alternatively, devTcsv utility allows for constructing that list from csv files

#these table have dependencies (cross references). However some may be compatible with others if cross references are. This allows for recombination of parameters

# A.1 axeT : description of  plant variability within the canopy
#-----------------------------------------------------------------------
#
# This table characterises the plant to plant variability within the canopy. If only one plant is given, then adel will clone that plant. For wheat, to have a correct simulation of tiller dynamics, a minimum of 30 plants is recommended.
# The table contents explicit parameters (eg number of leaves on axes) and indexes that reference a particular developmental pattern (parameter set) described in the other tables (dimT and phenT) to allows factorisation.
#
# There is one line per axe (see axeTdemo.csv). Columns are :
#
# "plant" : the plant number the axe belongs to
# "axe" : axe number (starting at 0 for main stem) within the plant
# "nf" : final number of leaves on the axe
# "end" : NA or date of end of growth of the axe for regressing ones
# "disp" : NA or date of disparition the axe
# "dimIndex" : reference (int) to a parameter set in dimT 
# "phenIndex" : reference to a parameter set in phenT
# "earIndex" : reference to a parameter set in earT, or NA if the axe does not bear an ear
# "emf1","ligf1","senf1","dispf1" : date of tip emergence, collar emergence, end of senescence and disparition of the first leaf of the axe. This allows for parameterising heterogeneous emergence date
#
#A.2 dimT : description of dimensions of axes (or a subset of axes)
#------------------------------------------------------------------
#
#Dimensions could be given for all axes (using one dimIndex par axe) or a subset of axes. Dimension are given for relative position within the axe (=position/final number of leaves). Actual dimension of axe are hence dependent on dimIndex and nf.
#
#TO BE DONE : add scale factor (or max dim) in axeT and describe here double normalised functions)
#
#There is one line per phytomer (see dimTdemo.csv) of all indexed axes. Columns are :
#
# "index" : the index refered to in axeT
# "nrel" : normalised phytomer position, starting from the base)
# "Ll": blade length
# "Lw": blade max width
# "Gl" : sheath length
# "Gd" : sheath diameter
# "El" : internode length
# "Ed" : internode diameter
#
#
#A.3 phenT : description of phenology of axes (or a subset of axes)
#------------------------------------------------------------------
#
# Phenology controls the rate of plant development (hence extension rates of organs), the dynamics of leaf appearance and the dynamics of senescence. 
#
#positions are normalised to final leaf numbers to allows sharing of data between axes of axeT table.
#
#dates of developmental events are given relative taking as origin the date of the event on leaf 1 of the axe. Actual development is computed from this table and the date concerning leaf 1 in axeT. 
#
#There is one line per phytomer (see phenTdemo.csv) of all indexed axes. Columns are :
#
# "index" : the index referred to in axeT
# "nrel" : normalised phytomer position, starting from 0 (to allow extrapolation)
# "tip" : date (origin tip leaf 1) of tip emergence of the phytomer
# "col" : date (origin col leaf 1) of collar emergence of the phytomer
# "ssi" : date (origin sen leaf 1) of full senescence of the phytomer (ssi)
# "disp" : date (origin disp leaf 1) of leaf disappearance.Blade disappear at disp. Sheath disappear when leaf above it disappear
#
# A.4 earT : description of ear dimension and phenology
# -----------------------------------------------------
#
#There is one line per ear type (refered by ear Index in axeT)
#
# "index" : the index refered to in axeT
# "em_ear" : delay between flag leaf ligulation and ear (tip of highest spike without awn) appearance
# "em_ped" : delay between flag leaf ligulation and peduncle (tip = base of th ear) appearance
# "end_gf" : delay between flag leaf ligulation and end of grain filling (full senescence of the ear+stem)
# "l_ped" : length of the peduncle
# "d_ped" : diameter of the peduncle
# "l_ear" : length of the ear (without awns)
# "Sp_ear" : projected area of ear without awn
# "l_ear_awn" : length of the ear+awns
#
#
#
# A 5: ssi2sen : desciption of progression of upper leaf senescence as a function of ssi
#
#by default, leaves start senescence 1 ssi unit before ssi = leaf number and comple senecence when ssi = leaf number.
#For the 'ndel' upper leaves, senecence start at sssi = final leaf number - ndel (t0) at a slower rate specified in this file (rate), and accelerate dssit1 ssi unit after t0, and leaves are fully senesced dssit2 ssi unit after t0.
#
#the table allows for definition of rate, dssit1 and dssit2 for the ndel upper leaves
# ndel is given by the number of lines of the file
#
#
# B. Inputs describing geometry
# -------------------------------

#geometry of leaves is given two list of list of matrices describing midrib curvature and leaf width variation with distance to the base of the leaf.
#the first level in the list is for collection index
# the second level is for matrix index. see alea

#
#beside leaf shapes two list of R function should be provided as inputs.
#
#The first list should provide 3 R function of axe number (0 = main stem) that returns:
# for azT : the azimuth(deg) of the first leaf of the axe with reference to the azimuth of the parent leaf
# for incT : the inclination (deg) of the base of the tiller compared with main stem
# for dredT : the distance (at maturity) between tiller and main stem

#These functions could be created with 'genGeoAxe' node (with constraints) or freely defined with 'freeGeoAxe'. A sample definition may be :

geoAxe <- list(
               azT = function(a) {
                 ifelse(a == 0,
                        0,
                        75 + (runif(1) - .5) * 5)
               },
               incT = function(a) {
                 ifelse(a == 0,
                        runif(1) * 5,
                        82 + (runif(1) - .5) * 5)
               },
               dredT = function(a) {
                 ifelse(a == 0,
                        0,
                        runif(1) * 7)
               }
               )
#
#The second list should provide two Rfunctions with args =  axe number, leaf position and finalnumber of leaves on the axe (nf). Returned values should be :
#
# for azim : the azimuth (deg) of the leaf compared to the previous one
# for Lindex : the positional index (starting at 1) of the collection to use for leaf shape. In case of leaf shape given in the form of a python dict, Lindex is the nth (starting at 1) sorted keys of the dict
#
#This list could be generated by genGeoLeaf or freeGeoLeaf. A sample dfinition may be : 
geoLeaf <- list(
                Azim = function(a,n,nf) {0 * runif(1)},
                Lindex = function(a,n,nf) {nf - n + 1}
                )

# C. Inputs describing simulation
# -------------------------------

#time step is given as a list of date for which a mock-up is wished
# position of plants within the plot are given externally from adel to a planter.
