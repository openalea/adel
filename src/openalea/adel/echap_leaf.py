#
#       echap_leaf
#
#       Copyright 2006-2018 INRIA - CIRAD - INRA
#
#       File author(s): Christian Fournier <Christian.Fournier@supagro.inra.fr>
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       OpenAlea WebSite : http://openalea.gforge.inria.fr
#
###############################################################################

"""
A median leaf shape parameterised as a function of leaf age
Data from echap project (contact: Corinne.Robert@inra.fr, Christian.Fournier@inra.fr)
"""

import os
import numpy
import pandas

from openalea.adel.geometric_elements import Leaves

datadir = os.path.dirname(__file__) + "/echap_leaf_data/"


def sr_data(sr_model="Mercia"):
    def sr_reader(fn):
        header_row_srdb = ["rankclass", "s", "r"]
        return pandas.read_csv(
            fn, names=header_row_srdb, sep=",", index_col=False, skiprows=1, decimal="."
        )

    srdb = {
        "Mercia": datadir + "srdb_GrignonMercia2010.csv",
        "Tremie": datadir + "srdb_GrignonTremie2011.csv",
    }

    return sr_reader(srdb[sr_model])


def read_trajectories(fn):
    dat = pandas.read_csv(
        fn,
        names=["age", "x", "y", "lindex"],
        sep=",",
        index_col=False,
        skiprows=1,
        decimal=".",
    )
    ages = set(dat["age"])
    numage = numpy.array(sorted(list(ages)))
    agemed = (numage[1:] + numage[:-1]) / 2.0
    bins = [numage[0] - 1.0] + agemed.tolist() + [numage[-1] + 1.0]
    dat["age_class"] = pandas.cut(dat["age"], bins, labels=False)
    grouped = dat.groupby(["lindex", "age_class"])
    trajs = {
        k: [
            {
                a: grouped.get_group((k, a)).loc[:, ["x", "y"]]
                for a in set(dat["age_class"])
            }
        ]
        for k in set(dat["lindex"])
    }
    return trajs, bins


def median_leaf_trajectories(xy_model="Tremie_byleafclass"):
    """interpolated xy for upper/lower leaves every 0.5 HS on Grignon 2009-2010 data"""
    trajs = {
        "MerciaRht_byleafclass": "MerciaRht_Grignon2010_byleafclass",
        "Tremie_byleafclass": "Tremie_Boigneville2012_2013_byleafclass",
        "Soissons_byleafclass": "Soissons_Grignon2010_byleafclass",
        "MerciaRht": "MerciaRht_Grignon2010",
        "Tremie": "Tremie_Boigneville2012_2013",
        "Soissons": "Soissons_Grignon2010",
        "Tremie12": "Tremie12",
        "Tremie13": "Tremie13",
    }
    fn = {k: datadir + "median_leaf_trajectories_" + trajs[k] + ".csv" for k in trajs}
    return read_trajectories(fn[xy_model])


def geoLeaf(nlim=4, dazt=60, dazb=10, Lindex_base=1, Lindex_top=2):
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
              ntop = nf - n
              ifelse(ntop <= {ntoplim:d}, {top_class:d}, {base_class:d})
              }}
              )
    """
    return rcode.format(
        ntoplim=nlim,
        dazTop=dazt,
        dazBase=dazb,
        top_class=Lindex_top,
        base_class=Lindex_base,
    )


def echap_leaves(
    xy_model="Tremie_byleafclass",
    sr_model="Mercia",
    disc_level=7,
    top_leaves=4,
    dazt=60,
    dazb=10,
):
    gL = geoLeaf(nlim=top_leaves, dazt=dazt, dazb=dazb)
    trajs, bins = median_leaf_trajectories(xy_model)
    sr = sr_data(sr_model)
    sr["Lindex"] = sr["rankclass"]
    srdb = {k: v.loc[:, ["s", "r"]].to_dict("list") for k, v in sr.groupby("Lindex")}
    return Leaves(
        trajs, srdb, geoLeaf=gL, dynamic_bins=bins, discretisation_level=disc_level
    )
