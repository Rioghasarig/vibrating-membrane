"""Time integration of M q'' + C q' + K q = F(t)."""

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse.linalg import splu

from .fem import Membrane


@dataclass
class Result:
    """Output of ``simulate``.

    ``frames`` is sampled at ``t``, for animation. The other series have one
    entry per time step, at ``step_t``.
    """

    membrane: Membrane
    dt: float
    t: np.ndarray  # (n_frames,)
    frames: np.ndarray  # (n_frames, N)
    step_t: np.ndarray  # (n_steps + 1,)
    kinetic: np.ndarray
    potential: np.ndarray
    work_in: np.ndarray  # cumulative work done by the drives
    dissipated: np.ndarray  # cumulative energy removed by damping
    probes: list = field(default_factory=list)  # [(x, y), ...]
    probe_u: np.ndarray = None  # (n_steps + 1, n_probes)

    @property
    def energy(self):
        return self.kinetic + self.potential

    def max_slope(self):
        """Largest |grad u| in any saved frame. The theory needs this << 1."""
        m = self.membrane
        u = m.to_full(self.frames)[:, m.triangles]
        area2 = 2 * m.areas
        ux = np.sum(u * m._grad_x_times_2A, axis=2) / area2
        uy = np.sum(u * m._grad_y_times_2A, axis=2) / area2
        return float(np.sqrt(ux**2 + uy**2).max())


def simulate(m, q0=None, v0=None, drives=(), duration=0.03, dt=None, n_frames=300, probes=None):
    """Integrate the membrane from initial displacement q0 and velocity v0.

    Uses the Newmark average-acceleration method. It is unconditionally
    stable and adds no numerical damping. Non-uniform damping makes C not
    proportional to M, so the equations can't be decoupled into modes. The
    system matrix is constant, so it is factored once and each step is a
    single back-substitution.

    Parameters
    ----------
    m : Membrane
    q0, v0 : (N,) vectors, zero if omitted
    drives : sequence of F(t) functions, added together
    duration : simulated time [s]
    dt : time step [s], default half the time a wave takes to cross an element
    n_frames : number of displacement frames to keep
    probes : list of (x, y) points to record at every step
    """
    if dt is None:
        dt = 0.5 * m.h / m.wave_speed
    n_steps = max(1, int(np.ceil(duration / dt)))
    dt = duration / n_steps

    q = np.zeros(m.N) if q0 is None else np.array(q0, dtype=float)
    v = np.zeros(m.N) if v0 is None else np.array(v0, dtype=float)
    drives = list(drives)

    def total_force(t):
        F = np.zeros(m.N)
        for drive in drives:
            F += drive(t)
        return F

    M, K, C = m.M, m.K, m.C
    if probes is None:
        probes = [(0.37 * m.a, 0.21 * m.a)]
    probe_weights = np.array([m.point_weights(*p) for p in probes])

    solve = splu((M + 0.5 * dt * C + 0.25 * dt**2 * K).tocsc()).solve

    frame_steps = np.unique(np.linspace(0, n_steps, n_frames).round().astype(int))
    frames = np.zeros((len(frame_steps), m.N))
    frames[0] = q
    next_frame = 1

    kinetic = np.zeros(n_steps + 1)
    potential = np.zeros(n_steps + 1)
    work_in = np.zeros(n_steps + 1)
    dissipated = np.zeros(n_steps + 1)
    probe_u = np.zeros((n_steps + 1, len(probes)))

    F = total_force(0.0)
    Kq = K @ q
    a = splu(M.tocsc()).solve(F - C @ v - Kq)
    kinetic[0] = 0.5 * v @ (M @ v)
    potential[0] = 0.5 * q @ Kq
    probe_u[0] = probe_weights @ q

    for n in range(n_steps):
        F_next = total_force((n + 1) * dt)

        q_pred = q + dt * v + 0.25 * dt**2 * a
        v_pred = v + 0.5 * dt * a
        a_next = solve(F_next - C @ v_pred - K @ q_pred)
        q_next = q_pred + 0.25 * dt**2 * a_next
        v_next = v_pred + 0.5 * dt * a_next

        # These increments make the discrete energy balance hold exactly:
        # dE = dq . (F + F_next) / 2 - dt * v_mid . C v_mid
        v_mid = 0.5 * (v + v_next)
        work_in[n + 1] = work_in[n] + (q_next - q) @ (0.5 * (F + F_next))
        dissipated[n + 1] = dissipated[n] + dt * v_mid @ (C @ v_mid)

        q, v, a, F = q_next, v_next, a_next, F_next
        kinetic[n + 1] = 0.5 * v @ (M @ v)
        potential[n + 1] = 0.5 * q @ (K @ q)
        probe_u[n + 1] = probe_weights @ q

        if next_frame < len(frame_steps) and n + 1 == frame_steps[next_frame]:
            frames[next_frame] = q
            next_frame += 1

    return Result(
        membrane=m,
        dt=dt,
        t=frame_steps * dt,
        frames=frames,
        step_t=np.arange(n_steps + 1) * dt,
        kinetic=kinetic,
        potential=potential,
        work_in=work_in,
        dissipated=dissipated,
        probes=list(probes),
        probe_u=probe_u,
    )
