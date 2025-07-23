from time import sleep as _sleep
from openalea.mtg.algo import height, rank
import numpy as np

def sleep(obj, seconds=0.0):
    '''sleep(seconds)

    Delay execution for a given number of seconds.  The argument may be
    a floating point number for subsecond precision.'''
    
    _sleep(seconds)
    return obj,

def leaf_sectors(g, leaf_number=0, latence=350):
    """
    filter_leaf_sectors(mtg,target leaf number, latency)
    
    extract 2 list of sectors from the mtg : sectors of a target leaf (0 means all leaves whose age <= latency) and sectors whose age > latency
    """
    max_scale = g.max_scale()
    age = g.property('age')

    infectious = [l for l in g.vertices(scale=max_scale) if 'Leaf' in g.label(l) and age.get(l) >= latence]

    if leaf_number == 0:
        green = [l for l in g.vertices(scale=max_scale) if 'Leaf' in g.label(l) and 0< age.get(l,0) < latence]
    else:
        scale_metamer = max_scale-1
        tips = [v for v in g.vertices(scale=scale_metamer) if g.is_leaf(v)]
        green = []
        for vid in tips:
            for i in range(1,leaf_number):
                vid = g.parent(vid)
            green.extend(l for l in g.components(vid) if 'Leaf' in g.label(l) and 0 < age.get(l,0) < latence )

    return g, list(map(g.node,green)), list(map(g.node, infectious))

def leaf_sectors_by_number(g, target_leaf_number=1, source_leaf_number=4):
    max_scale = g.max_scale()

    infectious = []
    green = []
    scale_metamer = max_scale-1
    tips = [v for v in g.vertices(scale=scale_metamer) if g.is_leaf(v)]
    for vid in tips:
	rvid = vid
        for i in range(1,target_leaf_number):
            rvid = g.parent(rvid)
            green.extend(l for l in g.components(rvid) if 'Leaf' in g.label(l))
	rvid = vid
        for i in range(1,source_leaf_number):
            rvid = g.parent(rvid)
            infectious.extend(l for l in g.components(rvid) if 'Leaf' in g.label(l))

    return g, list(map(g.node,green)), list(map(g.node, infectious))

def compute_distance(target_sectors, source_sectors, distance_function):
    """ Apply function on each target sectors as a function of other sectors"""
    if (not target_sectors) and (not source_sectors):
        return None

    for n in target_sectors:
        distance_function(n, source_sectors)

def get_distances(g,filename):
    """
    extract and reshape in a dict distance computed
    """

    d = g.property('distance')
    vals ={}
    ages, distances, leaves, sectors=[],[],[],[]
    for k in d:
        if len(d[k][0])==0:
            continue
	age, distance = list(zip(*d[k]))
	ln = height(g, g.complex(k))+1
	sector=rank(g,k)+1
	nrow = len(age)
	ages.extend(age)
	distances.extend(distance)
	leaves.extend([ln]*nrow)
	sectors.extend([sector]*nrow)


    ages = np.array(ages)
    distances = np.array(distances)
    leaves = np.array(leaves)
    sectors = np.array(sectors)

    res =np.array((leaves, sectors, ages, distances*100)).T

    np.savetxt(filename, res, delimiter=';', fmt=['%d', '%d','%d','%d'])
    return filename, leaves, sectors, ages, distances

