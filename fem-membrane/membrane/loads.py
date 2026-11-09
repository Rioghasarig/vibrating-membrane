"""Initial conditions and driving forces.

Initial conditions return an (N,) vector of interior node values.

Drives return a function F(t) giving the (N,) generalized force vector of

    M q'' + C q' + K q = F(t),    F_i(t) = integral(p(x, y, t) phi_i dA)

The force comes from the virtual work of the pressure p on a displacement
du, which is integral(p du dA). membrane-fem.md stops at F = 0.
"""

import warnings

import numpy as np
from scipy.special import jn_zeros, jv


# Initial conditions


def strike(m, x0, y0, impulse=1e-3, radius=0.01):
    """Velocity of a mallet strike at (x0, y0).

    The impulse [N s] is spread over a Gaussian patch of the given radius
    [m]. It is normalized so the membrane's total momentum equals ``impulse``.
    An ideal point strike isn't representable on a mesh, so keep ``radius``
    at about two elements or more.
    """
    if radius < 1.5 * m.h:
        warnings.warn(
            f"strike radius {radius:g} m is under-resolved by element size "
            f"{m.h:g} m; increase n_rings or radius",
            stacklevel=2,
        )
    shape = m.interior_values(
        lambda x, y: np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * radius**2))
    )
    return impulse * shape / (m.sigma * (m.node_areas @ shape))


def gaussian_bump(m, x0, y0, height=1e-3, radius=0.02):
    """Displacement shaped like a smooth bump."""
    return m.interior_values(
        lambda x, y: height * np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * radius**2))
    )


def pluck(m, x0, y0, height=1e-3):
    """Cone-shaped displacement: pulled up at (x0, y0), straight to the rim."""
    r0 = np.hypot(x0, y0)
    if r0 >= m.a:
        raise ValueError("pluck point must be inside the membrane")

    def cone(x, y):
        dx, dy = x - x0, y - y0
        dist = np.hypot(dx, dy)
        # Distance from (x0, y0) to the rim along the ray through (x, y).
        with np.errstate(invalid="ignore", divide="ignore"):
            ux, uy = dx / dist, dy / dist
        along = x0 * ux + y0 * uy
        to_rim = -along + np.sqrt(along**2 + m.a**2 - r0**2)
        return np.where(dist > 0, height * (1 - dist / to_rim), height)

    return m.interior_values(cone)


def bessel_mode(m, mm, nn, amp=1e-3):
    """Exact mode (mm, nn): amp * J_mm(j_mm,nn r / a) cos(mm theta)."""
    k = jn_zeros(mm, nn)[-1] / m.a
    return m.interior_values(
        lambda x, y: amp * jv(mm, k * np.hypot(x, y)) * np.cos(mm * np.arctan2(y, x))
    )


# Waveforms: s(t) accepts a scalar or an array of times


def sine(freq, amp=1.0, phase=0.0):
    """Continuous tone."""
    return lambda t: amp * np.sin(2 * np.pi * freq * t + phase)


def tone_burst(freq, cycles=5, amp=1.0):
    """A Hann-windowed tone lasting ``cycles`` periods from t = 0."""
    length = cycles / freq

    def waveform(t):
        t = np.asarray(t, dtype=float)
        window = np.where(
            (t >= 0) & (t <= length), np.sin(np.pi * t / length) ** 2, 0.0
        )
        return amp * window * np.sin(2 * np.pi * freq * t)

    return waveform


def gaussian_pulse(t0, width, amp=1.0):
    """A single pulse centred at t0 with standard deviation ``width``."""
    return lambda t: amp * np.exp(-((np.asarray(t, dtype=float) - t0) ** 2) / (2 * width**2))


# Drives


def plane_wave(m, waveform, angle=0.0, azimuth=0.0, c_air=343.0):
    """Acoustic plane wave. ``waveform`` is the pressure [Pa] over time.

    ``angle`` is degrees away from the membrane normal, so 0 hits every point
    at once. ``azimuth`` is the in-plane travel direction in degrees. Pressure
    at (x, y) is s(t - d(x, y) / c_air), where d is the distance travelled
    along the membrane scaled by sin(angle).

    This is a prescribed one-sided load: the membrane doesn't radiate or
    scatter sound back into the air.
    """
    sin_angle = np.sin(np.radians(angle))
    azimuth = np.radians(azimuth)
    delay = sin_angle * (m.x * np.cos(azimuth) + m.y * np.sin(azimuth)) / c_air
    return lambda t: m.B_ia @ waveform(t - delay)


def point_force(m, x0, y0, waveform):
    """Concentrated force at (x0, y0). ``waveform`` is the force [N] over time."""
    weights = m.point_weights(x0, y0)
    return lambda t: waveform(t) * weights


def pressure_field(m, fn):
    """Arbitrary pressure fn(x, y, t) [Pa] at node coordinate arrays x, y."""

    def force(t):
        p = np.broadcast_to(np.asarray(fn(m.x, m.y, t), dtype=float), (m.n_nodes,))
        return m.B_ia @ p

    return force
