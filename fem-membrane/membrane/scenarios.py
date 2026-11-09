"""Named simulation setups for simulate.py.

To add one, write a function decorated with @scenario that takes the parsed
command-line options and returns a Scenario. Each option is None unless the
user set it, so use ``option(opts.x, default)``.
"""

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from . import loads
from .fem import Membrane

RADIUS = 0.15  # m
SCENARIOS = {}


@dataclass
class Scenario:
    membrane: Membrane
    title: str
    duration: float
    q0: np.ndarray = None
    v0: np.ndarray = None
    drives: list = field(default_factory=list)
    probes: list = None
    n_frames: int = 300


def scenario(name):
    def register(fn: Callable):
        SCENARIOS[name] = fn
        return fn

    return register


def option(value, default):
    return default if value is None else value


def _membrane(opts, damping):
    return Membrane(a=RADIUS, n_rings=opts.rings, damping=option(opts.damping, damping))


@scenario("strike")
def strike(opts):
    """Mallet strike away from the centre, with light damping."""
    m = _membrane(opts, damping=4.0)
    x0, y0 = option(opts.at, (0.5 * m.a, 0.2 * m.a))
    return Scenario(
        membrane=m,
        title="Strike",
        duration=option(opts.duration, 0.05),
        v0=loads.strike(m, x0, y0),
        probes=[(x0, y0), (-0.4 * m.a, 0.3 * m.a)],
    )


@scenario("mode")
def mode(opts):
    """One exact Bessel mode oscillating freely (--mn M N)."""
    m = _membrane(opts, damping=0.0)
    mm, nn = option(opts.mn, (1, 1))
    return Scenario(
        membrane=m,
        title=f"Mode ({mm},{nn})",
        duration=option(opts.duration, 3 / m.analytic_frequency(mm, nn)),
        q0=loads.bessel_mode(m, mm, nn),
        n_frames=180,
    )


@scenario("tone")
def tone(opts):
    """Continuous acoustic plane wave, at the fundamental unless --freq is set."""
    m = _membrane(opts, damping=4.0)
    freq = option(opts.freq, m.f01)
    wave = loads.sine(freq, amp=option(opts.pressure, 2.0))
    drive = loads.plane_wave(m, wave, option(opts.angle, 0.0), option(opts.azimuth, 0.0))
    return Scenario(
        membrane=m,
        title=f"Plane wave, {freq:.0f} Hz",
        duration=option(opts.duration, 0.1),
        drives=[drive],
        n_frames=400,
    )


@scenario("blast")
def blast(opts):
    """Short oblique pressure pulse sweeping across the membrane."""
    m = _membrane(opts, damping=4.0)
    pulse = loads.gaussian_pulse(0.002, 0.0003, amp=option(opts.pressure, 200.0))
    drive = loads.plane_wave(m, pulse, option(opts.angle, 60.0), option(opts.azimuth, 30.0))
    return Scenario(
        membrane=m,
        title="Pressure pulse",
        duration=option(opts.duration, 0.03),
        drives=[drive],
    )


@scenario("damped-patch")
def damped_patch(opts):
    """Strike with a damping pad on one side: non-uniform damping."""
    base = option(opts.damping, 4.0)
    pad_x, pad_radius = 0.5 * RADIUS, 0.25 * RADIUS

    def pattern(x, y):
        return base + 300.0 * (np.hypot(x - pad_x, y) < pad_radius)

    m = Membrane(a=RADIUS, n_rings=opts.rings, damping=pattern)
    x0, y0 = option(opts.at, (-0.4 * RADIUS, 0.0))
    return Scenario(
        membrane=m,
        title="Strike with damping pad",
        duration=option(opts.duration, 0.05),
        v0=loads.strike(m, x0, y0),
        probes=[(x0, y0), (pad_x, 0.0)],
    )
