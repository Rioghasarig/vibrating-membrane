"""Checks of the simulator against results known from theory."""

import numpy as np
import pytest
from scipy.sparse.linalg import eigsh, spsolve
from scipy.special import jn_zeros, jv

from membrane import (
    Membrane,
    bessel_mode,
    gaussian_pulse,
    plane_wave,
    pluck,
    point_force,
    pressure_field,
    sine,
    simulate,
    strike,
    tone_burst,
)


def lowest_frequencies(m, k):
    eigenvalues = eigsh(m.K.tocsc(), k=k, M=m.M.tocsc(), sigma=0, which="LM")[0]
    return np.sqrt(np.sort(eigenvalues)) / (2 * np.pi)


# Mesh and matrices


def test_boundary_nodes_are_on_the_circle():
    m = Membrane(n_rings=12)
    assert np.allclose(np.hypot(m.x[m.boundary], m.y[m.boundary]), m.a)


def test_mesh_covers_the_inscribed_polygon():
    m = Membrane(n_rings=12)
    sides = len(m.boundary)
    polygon_area = 0.5 * sides * m.a**2 * np.sin(2 * np.pi / sides)
    assert m.areas.min() > 0
    assert m.areas.sum() == pytest.approx(polygon_area, rel=1e-12)


def test_matrices_are_symmetric_positive_definite():
    m = Membrane(n_rings=8)
    for matrix in (m.M, m.K):
        dense = matrix.toarray()
        assert np.allclose(dense, dense.T)
        assert np.linalg.eigvalsh(dense).min() > 0


# Spectrum


def test_frequencies_match_bessel_zeros():
    m = Membrane(n_rings=40)
    zeros = [
        jn_zeros(0, 1)[-1],
        jn_zeros(1, 1)[-1],
        jn_zeros(1, 1)[-1],
        jn_zeros(2, 1)[-1],
        jn_zeros(2, 1)[-1],
        jn_zeros(0, 2)[-1],
    ]
    exact = np.array(zeros) * m.wave_speed / (2 * np.pi * m.a)
    assert lowest_frequencies(m, 6) == pytest.approx(exact, rel=0.001)


def test_frequency_error_is_second_order_in_mesh_size():
    errors = []
    for rings in (10, 20, 40):
        m = Membrane(n_rings=rings)
        errors.append(abs(lowest_frequencies(m, 1)[0] - m.f01))
    assert errors[0] / errors[1] > 3
    assert errors[1] / errors[2] > 3


# Loads


def test_uniform_pressure_static_deflection():
    m = Membrane(n_rings=40)
    P = 10.0
    q = spsolve(m.K.tocsc(), pressure_field(m, lambda x, y, t: P)(0.0))
    r = np.hypot(m.x[m.interior], m.y[m.interior])
    exact = P * (m.a**2 - r**2) / (4 * m.T)
    assert np.abs(q - exact).max() < 0.01 * exact.max()


def test_normal_plane_wave_is_uniform_pressure():
    m = Membrane(n_rings=10)
    wave = plane_wave(m, sine(100.0, amp=3.0))
    uniform = pressure_field(m, lambda x, y, t: 3.0 * np.sin(2 * np.pi * 100.0 * t))
    assert wave(0.0013) == pytest.approx(uniform(0.0013))


def test_oblique_plane_wave_is_delayed_across_the_membrane():
    m = Membrane(n_rings=10)
    pulse = gaussian_pulse(0.0, 1e-4)
    F = plane_wave(m, pulse, angle=90.0, azimuth=0.0, c_air=100.0)
    # At t < 0 the wave has reached the -x side only. At t > 0, the +x side.
    x = m.x[m.interior]
    early = F(-m.a * 0.5 / 100.0)
    late = F(m.a * 0.5 / 100.0)
    assert early[x < 0].sum() > 100 * early[x > 0].sum()
    assert late[x > 0].sum() > 100 * late[x < 0].sum()


def test_point_force_total_matches_force():
    m = Membrane(n_rings=10)
    F = point_force(m, 0.03, 0.01, lambda t: 2.0)(0.0)
    assert F.sum() == pytest.approx(2.0)  # the point is away from the boundary


def test_point_outside_the_membrane_is_rejected():
    with pytest.raises(ValueError):
        Membrane(n_rings=5).point_weights(1.0, 0.0)


def test_strike_delivers_the_requested_impulse():
    m = Membrane(n_rings=30)
    v0 = strike(m, 0.02, 0.0, impulse=0.003, radius=0.01)
    assert m.sigma * m.node_areas @ v0 == pytest.approx(0.003)


def test_pluck_peaks_at_the_pluck_point_and_stays_positive():
    m = Membrane(n_rings=20)
    q = pluck(m, 0.0, 0.0, height=1e-3)
    assert q.max() == pytest.approx(1e-3)
    assert q.min() >= 0


def test_tone_burst_is_silent_outside_its_window():
    burst = tone_burst(100.0, cycles=4)
    assert burst(-0.001) == 0
    assert burst(0.041) == 0


# Time integration


def test_undamped_free_vibration_conserves_energy():
    m = Membrane(n_rings=15)
    result = simulate(m, v0=strike(m, 0.03, 0.02, radius=0.02), duration=0.02)
    assert np.ptp(result.energy) / result.energy[0] < 1e-9


def test_mode_oscillates_at_the_exact_frequency():
    m = Membrane(n_rings=30)
    amp = 1e-3
    q0 = bessel_mode(m, 0, 1, amp=amp)
    duration = 3 / m.f01
    result = simulate(m, q0=q0, duration=duration, probes=[(0.0, 0.0)])
    exact = amp * np.cos(2 * np.pi * m.f01 * result.step_t)
    assert np.abs(result.probe_u[:, 0] - exact).max() < 0.03 * amp


def test_damping_matrix_is_consistent():
    # Uniform c: C = integral(c phi_i phi_j) = (c / sigma) M exactly.
    m = Membrane(n_rings=8, damping=5.0)
    assert np.allclose(m.C.toarray(), (5.0 / m.sigma) * m.M.toarray())
    # C couples neighbouring nodes: it is not diagonal.
    assert m.C.nnz == m.M.nnz
    # No damping gives an empty matrix.
    assert Membrane(n_rings=8).C.nnz == 0


def test_uniform_damping_decays_at_the_expected_rate():
    c = 10.0
    m = Membrane(n_rings=30, damping=c)
    q0 = bessel_mode(m, 0, 1)
    duration = 0.1
    result = simulate(m, q0=q0, duration=duration)
    window = int(1 / (m.f01 * result.dt))  # average out the oscillation at 2 omega
    assert result.energy[-window:].mean() / result.energy[:window].mean() == pytest.approx(
        np.exp(-c * (duration - window * result.dt) / m.sigma), rel=0.02
    )


def test_energy_balance_with_driving_and_nonuniform_damping():
    m = Membrane(n_rings=15, damping=lambda x, y: 20.0 * (x > 0))
    wave = plane_wave(m, tone_burst(300.0, cycles=3, amp=5.0), angle=40.0)
    result = simulate(m, v0=strike(m, -0.03, 0.0, radius=0.02), drives=[wave], duration=0.03)
    balance = result.energy[0] + result.work_in - result.dissipated
    assert result.dissipated[-1] > 0
    assert result.work_in[-1] != 0
    assert np.abs(result.energy - balance).max() < 1e-9 * result.energy.max()


def test_steady_state_matches_the_closed_form_solution():
    """Uniform harmonic pressure on a damped membrane, measured at the centre.

    Z(r) = (P / D) (1 - J0(kappa r) / J0(kappa a)), where D = -sigma w^2 + i c w
    and kappa^2 = -D / T. This satisfies D Z - T lap(Z) = P with Z(a) = 0.
    """
    c, P, freq = 40.0, 2.0, 150.0
    m = Membrane(n_rings=24, damping=c)
    w = 2 * np.pi * freq
    drive = pressure_field(m, lambda x, y, t: P * np.cos(w * t))

    # Transients decay as exp(-c t / (2 sigma)), so this leaves about 1e-5.
    result = simulate(m, drives=[drive], duration=0.15, probes=[(0.0, 0.0)])

    D = -m.sigma * w**2 + 1j * c * w
    kappa = np.sqrt(-D / m.T + 0j)
    Z = (P / D) * (1 - 1 / jv(0, kappa * m.a))

    last_period = result.step_t > result.step_t[-1] - 1 / freq
    t = result.step_t[last_period]
    exact = np.real(Z * np.exp(1j * w * t))
    assert np.abs(result.probe_u[last_period, 0] - exact).max() < 0.02 * abs(Z)


def test_mode_shape_stays_fixed_at_the_rim_and_frames_are_recorded():
    m = Membrane(n_rings=12)
    result = simulate(m, q0=bessel_mode(m, 2, 1), duration=0.005, n_frames=11)
    assert len(result.t) == len(result.frames) == 11
    assert result.t[0] == 0
    assert result.t[-1] == pytest.approx(0.005)
    assert np.all(m.to_full(result.frames)[:, m.boundary] == 0)
