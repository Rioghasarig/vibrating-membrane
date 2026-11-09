"""Mass, stiffness, load and damping matrices (membrane-fem.md, sections 1-4)."""

import numpy as np
import scipy.sparse as sp
from scipy.special import jn_zeros

from .mesh import disk_mesh


class Membrane:
    """A fixed circular membrane discretized with linear triangles.

    Parameters
    ----------
    a : radius [m]
    sigma : surface mass density [kg/m^2]
    T : tension [N/m]
    n_rings : mesh resolution; element size is about a / n_rings
    damping : per-area damping coefficient c [N s/m^3]. None for no damping,
        a number for uniform damping, or a function c(x, y) of coordinate
        arrays for a non-uniform pattern. A function is evaluated at the
        triangle centroids.

    The damping matrix is C_ij = integral(c phi_i phi_j), from the virtual
    work of the damping force -c u_t (solver-equations.md, section 3.8). c is
    taken as constant on each triangle, so uniform c gives C = (c / sigma) M
    exactly, and every mode's amplitude decays as exp(-c t / (2 sigma)).

    Attributes
    ----------
    M, K, C : sparse (N, N) mass, stiffness and damping matrices over the
        N interior nodes
    B_ia : sparse (N, n_nodes) matrix turning nodal pressures into forces
    """

    def __init__(self, a=0.15, sigma=0.26, T=2000.0, n_rings=40, damping=None):
        self.a = a
        self.sigma = sigma
        self.T = T
        self.n_rings = n_rings
        self.nodes, self.triangles, self.boundary = disk_mesh(a, n_rings)
        self.x, self.y = self.nodes.T
        self.n_nodes = len(self.nodes)

        is_interior = np.ones(self.n_nodes, dtype=bool)
        is_interior[self.boundary] = False
        self.interior = np.flatnonzero(is_interior)
        self.N = len(self.interior)

        self._assemble(damping)

    def _assemble(self, damping):
        x = self.x[self.triangles]
        y = self.y[self.triangles]

        # grad(phi_i) = (b_i, c_i) / (2 A) on each triangle, for i = 0, 1, 2.
        b = np.column_stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]])
        c = np.column_stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]])
        self._grad_x_times_2A, self._grad_y_times_2A = b, c
        self.areas = 0.5 * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])

        # K^e_ij = T * A * grad(phi_i) . grad(phi_j)
        stiffness = (
            self.T
            * (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :])
            / (4 * self.areas[:, None, None])
        )
        # M^e_ij / sigma = integral(phi_i phi_j) = A (1 + delta_ij) / 12
        unit_mass = (self.areas[:, None, None] / 12) * (np.ones((3, 3)) + np.eye(3))

        rows = np.repeat(self.triangles, 3, axis=1).ravel()
        cols = np.tile(self.triangles, (1, 3)).ravel()
        shape = (self.n_nodes, self.n_nodes)
        K_all = sp.coo_matrix((stiffness.ravel(), (rows, cols)), shape=shape).tocsr()
        B_all = sp.coo_matrix((unit_mass.ravel(), (rows, cols)), shape=shape).tocsr()

        # Dropping boundary rows and columns gives the N x N matrices of the note.
        self.K = K_all[self.interior][:, self.interior]
        self.M = self.sigma * B_all[self.interior][:, self.interior]

        # F_i = integral(p phi_i). Rows are interior test functions, but the
        # columns keep every node because p need not vanish on the boundary.
        self.B_ia = B_all[self.interior]
        self.node_areas = np.asarray(self.B_ia.sum(axis=1)).ravel()

        # C^e_ij = integral(c phi_i phi_j) = c_e A (1 + delta_ij) / 12, with c
        # taken as its value c_e at the triangle's centroid.
        c_elem = self._element_damping(damping, x.mean(axis=1), y.mean(axis=1))
        damping_values = c_elem[:, None, None] * unit_mass
        C_all = sp.coo_matrix((damping_values.ravel(), (rows, cols)), shape=shape).tocsr()
        self.C = C_all[self.interior][:, self.interior]

    @staticmethod
    def _element_damping(damping, xc, yc):
        """Damping coefficient c on each triangle, from its centroid (xc, yc)."""
        if damping is None:
            return np.zeros(len(xc))
        c = damping(xc, yc) if callable(damping) else damping
        return np.broadcast_to(np.asarray(c, dtype=float), xc.shape)

    @property
    def wave_speed(self):
        return np.sqrt(self.T / self.sigma)

    @property
    def h(self):
        """Characteristic element size."""
        return self.a / self.n_rings

    def analytic_frequency(self, m, n):
        """Frequency [Hz] of the exact (m, n) mode: j_mn c / (2 pi a)."""
        return jn_zeros(m, n)[-1] * self.wave_speed / (2 * np.pi * self.a)

    @property
    def f01(self):
        """Exact fundamental frequency [Hz]."""
        return self.analytic_frequency(0, 1)

    def interior_values(self, fn):
        """Evaluate fn(x, y) at the interior nodes."""
        values = fn(self.x[self.interior], self.y[self.interior])
        return np.broadcast_to(np.asarray(values, dtype=float), (self.N,)).copy()

    def to_full(self, q):
        """Insert the zero boundary values. ``q`` is (..., N), result (..., n_nodes)."""
        q = np.asarray(q)
        full = np.zeros(q.shape[:-1] + (self.n_nodes,))
        full[..., self.interior] = q
        return full

    def point_weights(self, x0, y0):
        """Return w such that w . q is the displacement at (x0, y0).

        The same vector is the load produced by a unit point force there,
        because both are the hat functions evaluated at the point.
        """
        x = self.x[self.triangles]
        y = self.y[self.triangles]
        area2 = 2 * self.areas
        # Barycentric coordinates, which are the hat functions on each triangle.
        bary = np.column_stack(
            [
                self._grad_x_times_2A[:, i] * (x0 - x[:, (i + 1) % 3])
                + self._grad_y_times_2A[:, i] * (y0 - y[:, (i + 1) % 3])
                for i in range(3)
            ]
        ) / area2[:, None]

        best = np.argmax(bary.min(axis=1))
        if bary[best].min() < -1e-9:
            raise ValueError(f"point ({x0}, {y0}) is outside the mesh")

        full = np.zeros(self.n_nodes)
        full[self.triangles[best]] = bary[best]
        return full[self.interior]
