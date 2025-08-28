def CanMTGPlanter(g, positions, random_seed=0, azimuths=None):
    """arrange plant in a stand"""
    if hasattr(g, "planter"):
        g.planter(positions, random_seed, azimuths)
    else:
        from openalea.adel.mtg import planter

        planter(g, positions, random_seed)
    return (g,)
