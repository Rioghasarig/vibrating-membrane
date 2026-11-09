"""Finite element simulation of a fixed circular membrane.

Solves the equations derived in ../membrane-fem.md, plus a forcing term:

    M q'' + C q' + K q = F(t)

See ../solver-equations.md for the derivation.
"""

from .fem import Membrane
from .loads import (
    bessel_mode,
    gaussian_bump,
    gaussian_pulse,
    plane_wave,
    pluck,
    point_force,
    pressure_field,
    sine,
    strike,
    tone_burst,
)
from .solver import Result, simulate

__all__ = [
    "Membrane",
    "Result",
    "simulate",
    "strike",
    "gaussian_bump",
    "pluck",
    "bessel_mode",
    "sine",
    "tone_burst",
    "gaussian_pulse",
    "plane_wave",
    "point_force",
    "pressure_field",
]
