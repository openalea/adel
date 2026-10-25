"""
The shoot's colours for metafspm's matplotlib drawings (openalea.metafspm.display.draw): its elements green, or straw
when more senesced than green.

    from openalea.metafspm.display import draw
    draw(axes, architecture.triangles(), colors=state_colors(ds.get("is_green")))
"""
import numpy as np

GREEN, SENESCED = (0.2, 0.55, 0.15), (0.8, 0.7, 0.35)


def state_colors(is_green):
    """Per element: green, or straw when it is more senesced than green."""
    return np.where(np.asarray(is_green, dtype=bool)[:, None], GREEN, SENESCED)
