import os
import pandas


#
# create AdelR functions into rpy2 environment
import numpy
import rpy2.robjects as robj
from rpy2.robjects import numpy2ri


r = robj.r

try:
    robj.globalEnv = robj.globalenv
except:
    pass


dir = os.path.dirname(__file__)

# content of rfiles using setuptools pkg_resources utility
# from pkg_resources import resource_string

# Set Numeric Locale value
r('Sys.setlocale(category="LC_NUMERIC",locale="C")')


def get_rcode(file_name):
    with open(os.path.join(dir, file_name), "r") as content_file:
        content = content_file.read()
    return content


# r.source(os.path.join(path,'Adel.R')) : code qui suit plutot a laisser au niveau des fonctions
rcode = get_rcode("Adel.R")
r(rcode)
rcode = get_rcode("setAdel.R")
r(rcode)
rcode = get_rcode("UseAdel.R")
r(rcode)
rcode = get_rcode("ArvalisToAdel.R")
r(rcode)
rcode = get_rcode("genString.R")
r(rcode)
# r.source(os.path.join(path,'ArvalisToAdel.R'))
RrunAdel = robj.globalEnv["runAdel"]
RsetCanopy = robj.globalEnv["setCanopy"]
RsetAdel = robj.globalEnv["setAdeluser"]
RdevCsv = robj.globalEnv["devTcsv"]
RreadCsv = robj.globalEnv["readCsv"]
# RgenGeoAxe = robj.globalEnv['genGeoAxe']
# RgenGeoLeaf = robj.globalEnv['genGeoLeaf']
RsetAdelArv = robj.globalEnv["setAdelArv"]
RgenString = robj.globalEnv["genString"]
RcanL2canS = robj.globalEnv["canL2canS"]
RcheckAxeDyn = robj.globalEnv["checkAxeDyn"]
RgetAxeT = robj.globalEnv["getAxeT"]
RgetPhenT = robj.globalEnv["getPhenT"]
RgetPhytoT = robj.globalEnv["getPhytoT"]
# RgetLeafT = robj.globalEnv['getLeafT']


def readRData(fn):
    """return a dictionary containing the Robject of the RData file"""

    data = r.load(fn)
    return dict(list(zip([nom for nom in data], [r[nom] for nom in data])))


def saveRData(Robj, name, fn):
    """name and write the Robj to RData file"""
    r.assign(name, Robj)
    r.save(list=name, file=fn)
    return fn


def RlistAsDict(Rlist):
    """returns a dictionary containing the elements of the Rlist"""
    if r["is.null"](Rlist.names)[0]:
        Rlist.names = r["seq"](Rlist)
    return dict(list(zip([n for n in r.names(Rlist)], [obj for obj in Rlist])))


def _rvect_asarray(rvect):
    """ "
    convert an r_vector into an array of numeric values or into an array of string if NA are present (numpy2ri replaces NA with numeric values otherwise) .
    Will be deprecated when numpy will offer true NA type"""

    if r["length"](r["which"](r["is.na"](rvect)))[0] > 0:
        return numpy.array(
            r["ifelse"](r["is.na"](rvect), "NA", rvect)
        )  # as.character fails if first value is less than one character  (eg (1,NA,..))
    elif isinstance(rvect, robj.vectors.FactorVector):
        return numpy.array(r["as.character"](rvect))
    else:
        return numpy.array(rvect)


def dataframeAsdict(df):
    """convert an RDataframe to a python dict"""
    if r["is.null"](df)[0]:
        return None
    #    try:
    #        d = dict(zip( df.colnames, numpy.array(df)))
    #        return d
    #    except:
    if isinstance(df, numpy.recarray):
        d = {k: _rvect_asarray(df[k]) for k in df.dtype.names}
    else:
        try:
            d = dict([(k, _rvect_asarray(df.r[k][0])) for k in r.colnames(df)])
        except:
            d = dict(
                [(k, _rvect_asarray(df.rx2(str(k)))) for k in r.colnames(df)]
            )  # r delegator is replaced by rx/rx2 in new rpy2
    return d


def _is_iterable(x):
    try:
        x = iter(x)
    except TypeError:
        return False
    return True


def dataframe(d):
    """convert a dict of numbers to an RDataframe"""
    cv = numpy2ri.converter
    with cv.context():
        df = {}
        if d is None:
            return r("as.null()")
        else:
            for k, v in d.items():
                rval = numpy.array(v)
                if not _is_iterable(v):
                    v = [v]
                if "NA" in numpy.array(v).tolist():
                    df[k] = r["as.numeric"](rval)
                else:
                    df[k] = rval
        dataf = r["data.frame"](**df)
        return dataf


def Rdflist(dictOfdict):
    """convert a dict of dict of numpy vectors into a Rlist of Rdataframe"""
    cv = numpy2ri.converter + robj.default_converter
    with cv.context():
        df_tags = list(dictOfdict.keys())
        df_value = list(dictOfdict.values())
        return r.list(**dict(list(zip(df_tags, [dataframe(dft) for dft in df_value]))))


def RdflistAsdicts(Rdflist):
    """convert a Rlist of Rdataframe into a dict of dict of numpy vectors"""
    return dict(
        list(
            zip(
                [n for n in r.names(Rdflist)], [dataframeAsdict(obj) for obj in Rdflist]
            )
        )
    )


def csvAsDict(fn, type=1):
    """returns a dictionnary with the content of csv file as numpy vectors (one colum = one key)"""
    df = RreadCsv(fn, type)
    return dataframeAsdict(df)


def setAdelArv(Rcalage, Rfunstr, np, sdlev=20):
    """Creates a set of parameter for simulating np plants with adel from arvalis calage Rdata (deprecated)"""

    r(Rfunstr)
    RFun = robj.globalEnv["Rfun"]
    p = RsetAdelArv(Rcalage, np, sdlev, RFun)
    return p


def setAdel(
    devT,
    RcodegeoLeaf,
    RcodegeoAxe,
    nplants=1,
    seed=None,
    xydb=None,
    srdb=None,
    sample="random",
    ssipars=None,
):
    """Creates a set of parameter for simulating np plants with adel from R inputs (see adeldoc.R)"""
    if seed is None:
        rseed = r("as.null()")
    else:
        rseed = seed

    if xydb is None:
        rxydb = r("as.null()")
    elif isinstance(
        xydb, dict
    ):  # pure python xydb dict, a list of list is passed to setAdel
        # Warning : setadel uses index only and finding the closest index if db is too small.
        # Unambiguous meaning of keys retrieving from Lindex afterward therefore will need sorting keys
        # before using Lindex (python sorted key index = Lindex - 1 (Lindex follows Rindexing convention),
        # python list index = lseed - 1. cf conversion infra
        keys = list(xydb.keys())
        keys.sort()  # to ensure lseed index is sample in the good list
        rxydb = r.list(*list(map(str, keys)))
        for i in range(len(rxydb)):
            rxydb[i] = r.list(*list(range(len(xydb[keys[i]]))))
    else:
        rxydb = r.load(xydb)[0]
        rxydb = r(rxydb)

    if srdb is None:
        rsrdb = r("as.null()")
    elif isinstance(srdb, dict):  # pure python srdb dict, a list is passed to setAdel
        rsrdb = r.list(*list(map(str, list(srdb.keys()))))
    else:  # path to RData file
        rsrdb = r.load(srdb)[0]
        rsrdb = r(rsrdb)

    RdevT = Rdflist(devT)
    r(RcodegeoAxe)
    geoAxe = robj.globalEnv["geoAxe"]
    r(RcodegeoLeaf)
    geoLeaf = robj.globalEnv["geoLeaf"]

    if ssipars is None:
        ssipars = r("as.null()")
    else:
        ssipars = robj.r["list"](**ssipars)
    p = RsetAdel(RdevT, geoLeaf, geoAxe, nplants, sample, rseed, rxydb, rsrdb, ssipars)
    return p


def leaf_keys(lindex, lseed, db):
    """convert R-style lindex/lseed (also called LcType/Lindex in canopy table)
    into (keys,index) of python xy/sr databases
    """
    if 1 > lindex or lindex > len(db) or lseed < 1:
        raise KeyError("invalid index for leaf shape database")
    keys = list(db.keys())
    keys.sort()
    return keys[lindex - 1], lseed - 1


def plantSample(setAdelPars):
    """return id of plants used by setAdel to set up the canopy"""
    p = r.names(setAdelPars)
    p = r["as.numeric"](p)
    return _rvect_asarray(p)


def checkAxeDyn(setAdelPars, dates, plant_density=1):
    """Compute the number of axes present in the canopy at given dates"""
    if type(dates) is not list:
        dates = [dates]
    d = robj.FloatVector(dates)
    df = RcheckAxeDyn(d, setAdelPars, plant_density)
    return pandas.DataFrame(dataframeAsdict(df))


def getAxeT(setAdelPars):
    """return axeT table"""
    df = RgetAxeT(setAdelPars)
    return pandas.DataFrame(dataframeAsdict(df))


def getPhenT(setAdelPars, axe="MS"):
    """return phenT table"""
    df = RgetPhenT(setAdelPars, axe=axe)
    return pandas.DataFrame(dataframeAsdict(df))


def getPhytoT(setAdelPars, axe="MS"):
    """return phytoT table"""
    df = RgetPhytoT(setAdelPars, axe=axe)
    return pandas.DataFrame(dataframeAsdict(df))


# def getLeafT(setAdelPars):
#    """ return a dict of leaves plantn_axetype_metamer : argument_dict_ for_leaf_sampler"""
#    def id_args(label):
#        plant, axe, n, nf = label.split('_')
#        id = '_'.join((plant, axe, n))
#        args = {'axe':axe, 'n': int(n), 'nf': int(nf)}
#        return (id,args)
#    return dict([id_args(s) for s in RgetLeafT(setAdelPars)])


def canL2canS(RcanT, srdb, shrink=1):
    """outputs are :
    - Ll : length of the blade
    - Lv : visible (emerged) length of blade (green + senesced, rolled + unrolled)
    - Lr : rolled part of the blade
    - Lsen : length of the senescent part of the blade (hidden + visible)
    - L_shape : Mature length of the blade used to compute blade shape
    - Lw_shape : Maximal width of the blade used to compute blade shape
    - Linc : relative inclination of the base of leaf blade at the top of leaf sheath (deg)
    - Laz : relative azimuth of the leaf
    - Lsect : the number of sectors per leaf blade
    - Gl : length of the sheath (hidden + visible)
    - Gv : emerged length of the sheath
    - Gsen : senescent length of the sheath (hidden + visible)
    - Gd : apparent diameter of the sheath
    - Ginc : relative inclination of the sheath
    - El: length of the internode (hidden + visible)
    - Ev: emerged length of the internode
    - Esen: senescent length of the internode (hidden + visible)
    - Ed: diameter of the internode
    - Einc : relative inclination of the internode

    """

    sr = r.load(srdb)[0]
    sr = r(sr)
    res = RcanL2canS(RcanT, sr, shrink)
    d = dataframeAsdict(res)
    return d


def setCanopy(RcanT, nplants=1, randomize=True, seed=None):
    """Duplicates or subsample a canopy table to create a new canopy with specified number of plants. Randomise allows for random sample and random azimuth"""
    if seed is None:
        rseed = r("as.null()")
    else:
        rseed = seed

    if randomize:
        rrand = 1
    else:
        rrand = 0

    can = RsetCanopy(RcanT, nplants, rrand, rseed)
    return can


def RunAdel(
    datesTT,
    plant_parameters,
        adelpars=None,
):
    """Run Adel model for each date in datesTT according to parameter list"""

    if adelpars is None:
        adelpars = {
            "senescence_leaf_shrink": 0.5,
            "leafDuration": 2,
            "fracLeaf": 0.2,
            "stemDuration": 2.0 / 1.2,
            "dHS_col": 0.2,
            "dHS_en": 0,
            "epsillon": 1e-6,
            "HSstart_inclination_tiller": 1,
            "rate_inclination_tiller": 30,
            "drop_empty": True,
        }
    if type(datesTT) is not list:
        datesTT = [datesTT]
    x = robj.FloatVector(datesTT)
    ap = robj.r["list"](**adelpars)
    # chn = RrunAdel(x,plant_parameters,ap)
    # return [c[0] for c in chn]
    res = RrunAdel(x, plant_parameters, ap)
    if len(res) <= 0:  # empty canopy
        d = None
    else:
        d = dataframeAsdict(res[0])
    return d


def _devCsv(axeTfn, dimTfn, phenTfn, earTfn=None, ssi2senTfn=None):
    """Import development parameters for adel from csv files and/or pandas dataframes"""
    args = [axeTfn, dimTfn, phenTfn]
    if earTfn is not None:
        args += [
            earTfn,
        ]
    if ssi2senTfn is not None:
        args += [
            ssi2senTfn,
        ]
    # use tempdir + csv to pass correctly NA from pandas to R, and to let RdevCsv make the columns renaming, if any
    if any([isinstance(f, pandas.DataFrame) for f in args]):
        import tempfile
        import shutil

        tmp_dir = tempfile.mkdtemp()
        for i in range(len(args)):
            df = args[i]
            if isinstance(df, pandas.DataFrame):
                f = tmp_dir + "/args%d.csv" % (i)
                df.to_csv(f, na_rep="NA", index=False)
                args[i] = f
        Rdat = RdevCsv(*args)
        shutil.rmtree(tmp_dir)
    else:
        Rdat = RdevCsv(*args)
    return Rdat


def devCsv(axeTfn, dimTfn, phenTfn, earTfn=None, ssi2senTfn=None):
    """Import development parameters for adel from csv files and/or pandas dataframes"""
    return RdflistAsdicts(
        _devCsv(axeTfn, dimTfn, phenTfn, earTfn=None, ssi2senTfn=None)
    )


def genString(RcanopyT):
    """Generate a Lsystem string from an R dataframe representing the canopy"""
    chn = RgenString(RcanopyT)
    return chn[0]


def genGeoAxe(
    azM=75, daz=5, ibmM=2, dibm=2, incT=60, dinT=5, dep=7, depMin=1.5, dazTb=60
):
    """generate geoAxe R code for Adel

    Parameters:
        - "axe Azimuth"(azM): define the rotation of a tiller around its axe, if azTb_true=0 and azM=90 then leaf1 of tiller is in the same plane as the beaning leaf
        - "dazT": is the random deviation around azM (azM_true = azM +/- daz/2)
        - "MainStem inclination"(incBm): is the  inclination of the main stem
        - "dincBm": is the random deviation around incBm (incBm_true = (incBm +/- dincBm/2)
        - "tiller inclination"(incT) is the basal inclination of the tiller relative to the main stem position
        - "dincT": is the random deviation around incT
        - "dredTMax": is the maximal distance at which the tiller will go upward (the distance at flowering between top of parent axe and top of tiller)
        - "dredMin": is the minimal distance at which tiller go upward
        - "dazTb": is the random deviation around "azTb" which is the azimuth angle between tiller axe and the midrib of bearing leaf. we suppose that (azTb = 0), so (azTb_true = 0 +/- dazTb/2)
    """
    rcode = """
    geoAxe <- list(
    azT = function(a) {{
          ifelse(a == 'MS',
          runif(1) * 360,#plant azimuth
          {azTM:.2f} + (runif(1) - .5) * {dazT:.2f})
       }},
    azTb = function(a) {{
         ifelse(a == 'MS',
                0,
                (runif(1) - .5) * {dazTb:.2f})
       }},
    incT = function(a) {{
       ifelse(a == 'MS',
              {incBmM:.2f} + (runif(1) - .5) * {dincBm:.2f},
              {incT:.2f} + (runif(1) - .5) * {dincT:.2f})
              }},
    dredT = function(a) {{
         #1.5 is an offset to avoid tiller superposed to mainstem
          ifelse(a == 'MS',
                 0,
                 {depMin:.2f} + runif(1) * ({depMax:.2f}-{depMin:.2f}))
        }}
       )
       """
    return rcode.format(
        azTM=azM,
        dazT=daz,
        incBmM=ibmM,
        dincBm=dibm,
        incT=incT,
        dincT=dinT,
        depMax=dep,
        depMin=depMin,
        dazTb=dazTb,
    )


def genGeoLeaf(nlim=4, dazt=60, dazb=10):
    """generate geoLeaf function for Adel"""
    rcode = """
    geoLeaf <- list(
     Azim = function(a,n,nf) {{
            ntop = nf - n
            ifelse(ntop <= {ntoplim:d},
            180 + {dazTop:.2f} * (runif(1) - .5),
            180 + {dazBase:.2f} * (runif(1) - .5))
            }},
     Lindex = function(a,n,nf) {{
              nf - n + 1}}
              )
        """
    return rcode.format(ntoplim=nlim, dazTop=dazt, dazBase=dazb)


def freeGeoAxe(rcode):
    """returns geoAxe code from txt definition"""

    return rcode


def freeGeoLeaf(rcode):
    """returns geoLeaf code from txt definition"""

    return rcode


def R_xydb(rdata):
    """return a python readable database from RData file containing leaf xy data
    - rdata is a regular .RData file containing a named list of list of xy matrix/dataframes
    """
    xy = r.load(rdata)[0]
    # rename Rlists with list index, to ensure leaf retrieval by index as key in python
    r("names(%s) = seq(%s)" % (xy, xy))

    xy = RlistAsDict(r(xy))

    rank = list(xy.keys())
    leaves = {}
    for k in rank:
        xyk = RlistAsDict(xy[k])
        leaves_id = list(xyk.keys())
        for leaf_id in leaves_id:
            try:
                x, y = numpy.transpose(numpy.array(xyk[leaf_id]))  # xy is a matrix
            except ValueError:  # xy is a dataframe
                x, y = numpy.array(xyk[leaf_id][0]), numpy.array(xyk[leaf_id][1])
            leaves.setdefault(k, []).append((x, y))
    return leaves


def R_srdb(rdata):
    """return a python readable database from RData file containing leaf sr data
    - rdata is a regular .RData file containing a named list of list of xy matrix/dataframes
    """
    sr = r.load(rdata)[0]
    # rename Rlists with list index, to ensure leaf retrieval by index as key in python
    r("names(%s) = seq(%s)" % (sr, sr))

    sr = RlistAsDict(r(sr))

    rank = list(sr.keys())

    leaves = {}
    for k in rank:
        try:
            s, radius = numpy.transpose(numpy.array(sr[k]))
        except ValueError:  # sr is a dataframe
            s, radius = numpy.array(sr[k][0]), numpy.array(sr[k][1])
        leaves[k] = (s, radius)

    return leaves
