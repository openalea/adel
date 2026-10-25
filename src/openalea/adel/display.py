"""
Matplotlib drawings of triangles (metafspm's coupling.geometry.Triangles, or their vertices), without PlantGL: a canopy
coloured by its elements' state (green or senesced) or by a value per element (light, ...).

    figure, axes = plt.subplots(subplot_kw=dict(projection="3d"))
    draw(axes, architecture.triangles(), colors=state_colors(ds.get("is_green")))
"""
import numpy as np

GREEN, SENESCED = (0.2, 0.55, 0.15), (0.8, 0.7, 0.35)


def state_colors(is_green):
    """Per element: green, or straw when it is more senesced than green."""
    return np.where(np.asarray(is_green, dtype=bool)[:, None], GREEN, SENESCED)


def value_colors(values, colormap="viridis", limits=None):
    """Per element: its value's colour in *colormap* (between *limits*, by default the values' range)."""
    from matplotlib import colormaps
    values = np.asarray(values, dtype=float)
    low, high = limits if limits is not None else (np.nanmin(values), np.nanmax(values))
    return colormaps[colormap](np.clip((values - low) / ((high - low) or 1.), 0., 1.))[:, :3]


def draw(axes, triangles, colors=None, owner=None, triangle_colors=None, view=(15., -60.), shade=True):
    """Draw *triangles* (Triangles, or vertices (T, 3, 3) with their *owner*) on 3d *axes*, each triangle coloured by
    its element's row of *colors* (per element, (n, 3)) or by its row of *triangle_colors* ((T, 3)); the axes scaled
    equally."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    vertices = np.asarray(getattr(triangles, "vertices", triangles), dtype=float)
    owner = getattr(triangles, "owner", owner)
    if triangle_colors is not None:
        facecolors = np.asarray(triangle_colors, dtype=float)
    elif colors is not None:
        facecolors = np.asarray(colors, dtype=float)[owner]
    else:
        facecolors = np.tile(GREEN, (len(vertices), 1))
    if shade and len(vertices):
        normal = np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0])
        norm = np.linalg.norm(normal, axis=1, keepdims=True)
        light = np.abs(normal[:, 2:3] / np.where(norm > 0, norm, 1.))         # lit from above
        facecolors = facecolors * (0.6 + 0.4 * light)
    axes.add_collection3d(Poly3DCollection(vertices, facecolors=facecolors, edgecolors="none"))
    if len(vertices):
        points = vertices.reshape(-1, 3)
        low, high = points.min(axis=0), points.max(axis=0)
        centre, half = (low + high) / 2., max((high - low).max() / 2., 1e-3)
        axes.set_xlim(centre[0] - half, centre[0] + half)
        axes.set_ylim(centre[1] - half, centre[1] + half)
        axes.set_zlim(max(0., centre[2] - half), max(0., centre[2] - half) + 2 * half)
    axes.set_box_aspect((1., 1., 1.))
    axes.view_init(*view)
    axes.set_xlabel("x (m)")
    axes.set_ylabel("y (m)")
    axes.set_zlabel("z (m)")
    return axes
