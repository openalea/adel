from openalea.mtg.io import lpy2mtg as _lpy2mtg


def lpy2mtg(axialtree, lsystem):
    """ """
    scene = lsystem.sceneInterpretation(axialtree)
    g = _lpy2mtg(axialtree, lsystem, scene=scene)
    geometry = g.property("geometry")
    geometry.clear()
    geoms = dict((s.id, s) for s in scene)
    vids = g.vertices(scale=g.max_scale())
    vids.sort()
    vids.reverse()

    gk = list(set(geoms))
    gk.sort()
    gk.reverse()

    gids = list(zip(vids, gk))
    for k, v in gids:
        geometry[k] = geoms[v]
    return (g,)
