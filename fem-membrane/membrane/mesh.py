"""Triangular mesh of a disk."""

import numpy as np
from scipy.spatial import Delaunay


def disk_mesh(a, n_rings):
    """Mesh the disk of radius ``a`` with linear triangles.

    Nodes sit on concentric rings: one centre node plus 6k nodes on ring k,
    which keeps the triangles close to equilateral. The outer ring lies
    exactly on the circle x^2 + y^2 = a^2.

    Returns
    -------
    nodes : (n_nodes, 2) float array
    triangles : (n_triangles, 3) int array, counter-clockwise
    boundary : (n_boundary,) int array of node indices on the circle
    """
    if n_rings < 1:
        raise ValueError("n_rings must be at least 1")

    rings = [np.zeros((1, 2))]
    for k in range(1, n_rings + 1):
        theta = 2 * np.pi * np.arange(6 * k) / (6 * k)
        r = a * k / n_rings
        rings.append(np.column_stack([r * np.cos(theta), r * np.sin(theta)]))
    nodes = np.vstack(rings)

    delaunay = Delaunay(nodes)
    if len(delaunay.coplanar):
        raise RuntimeError("Delaunay triangulation dropped mesh nodes")
    triangles = delaunay.simplices

    # Signed doubled area, used to remove slivers and fix orientation.
    p = nodes[triangles]
    area2 = (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1]) - (
        p[:, 2, 0] - p[:, 0, 0]
    ) * (p[:, 1, 1] - p[:, 0, 1])
    keep = np.abs(area2) > 1e-10 * (a / n_rings) ** 2
    triangles = triangles[keep]
    flip = area2[keep] < 0
    triangles[flip] = triangles[flip][:, [0, 2, 1]]

    boundary = np.arange(len(nodes) - 6 * n_rings, len(nodes))
    return nodes, triangles, boundary
