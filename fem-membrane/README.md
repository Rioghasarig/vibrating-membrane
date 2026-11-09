# FEM vibrating membrane

Simulates a fixed circular membrane with the finite element equations from [`../membrane-fem.md`](../membrane-fem.md):

    M q'' + C q' + K q = F(t)

[`solver-equations.md`](solver-equations.md) derives every term and lists the equations the code evaluates.

## Setup

Python 3.10 or newer. From this directory:

    pip install -r requirements.txt
    pytest

The tests compare the code against known theory, so run them first. `.gif` output works out of the box, and `.mp4` also needs ffmpeg on your `PATH`.

## Running

    python simulate.py strike
    python simulate.py mode --mn 2 1 --view 3d
    python simulate.py tone --angle 45 --freq 350
    python simulate.py blast --save blast.gif --save-summary blast.png --no-show

| Scenario | What it does |
|---|---|
| `strike` | Off-centre mallet strike with light damping. Move it with `--at X Y`. |
| `mode` | One exact Bessel mode oscillating freely. Choose it with `--mn M N`. |
| `tone` | Continuous acoustic plane wave, at the fundamental unless `--freq` is set. |
| `blast` | Short oblique pressure pulse sweeping across the membrane. |
| `damped-patch` | Strike with a damping pad on one side (non-uniform damping). |

Each run opens an animation and a summary figure with probe displacement, energy and spectrum. Use `--no-show` with `--save` or `--save-summary` if you have no display. `python simulate.py --help` lists every option.

The default membrane is drum-like: a = 0.15 m, sigma = 0.26 kg/m², T = 2000 N/m. That gives a wave speed of 88 m/s and a fundamental of 224 Hz. Animations play in slow motion, and the title shows the simulated time.

The `3d` view rebuilds the surface every frame, so it's slow live at high ring counts. Save it to a file, or lower `--rings`.

## Writing your own simulation

Combine any initial condition with any number of drives:

```python
import matplotlib.pyplot as plt
from membrane import Membrane, simulate, strike, plane_wave, tone_burst, point_force, sine
from membrane.viz import animate, plot_summary

m = Membrane(a=0.2, sigma=0.3, T=1500, n_rings=50, damping=3.0)

result = simulate(
    m,
    v0=strike(m, x0=0.05, y0=0.0),
    drives=[
        plane_wave(m, tone_burst(400, cycles=8, amp=5.0), angle=30),
        point_force(m, -0.05, 0.05, sine(250, amp=0.01)),
    ],
    duration=0.05,
    probes=[(0.05, 0.0), (0.0, 0.0)],
)

anim = animate(result, view="2d")  # keep this reference until it has played
plot_summary(result)
plt.show()
```

Initial conditions are `strike` (velocity, passed as `v0`), and `gaussian_bump`, `pluck` and `bessel_mode` (displacement, passed as `q0`).

Drives are `plane_wave`, `point_force` and `pressure_field`, which takes any `p(x, y, t)`. Waveforms for the first two are `sine`, `tone_burst`, `gaussian_pulse` or your own `s(t)`.

`damping` can also be a function `c(x, y)`, for example to damp one region only.

To add a named command-line scenario, write a short `@scenario("name")` function in [`membrane/scenarios.py`](membrane/scenarios.py).

## How it extends the note

**Forcing.** The note ends with no external force. A pressure `p` does virtual work `∫ p δu dA`, which gives the generalized force `F_i(t) = ∫ p φ_i dA`. The code interpolates `p` onto the mesh nodes, including boundary nodes, and multiplies by the unit mass matrix. A point force `f(t)` at `x0` gives `F_i = f(t) φ_i(x0)`.

**Damping.** This follows step 4 of the note. You give a damping force per unit area, `−c(x, y) u_t`, with `c` in N·s/m³. This gives the full damping matrix `C_ij = ∫ c φ_i φ_j dA`, which couples neighbouring nodes in the same way `M` does. `c` is taken as constant on each triangle, at its centroid value. For uniform `c`, `C = (c/σ) M` exactly and every mode decays as `exp(-c t / 2σ)`. `c` is the same quantity as `γ` in the modal notes.

## Numerical method

Time stepping is Newmark average acceleration. It's unconditionally stable and adds no numerical damping, but the step still has to resolve the motion. The default step is half the time a wave takes to cross one element. If you drive above a few kHz, pass a smaller `--dt` and raise `--rings`, since the mesh needs about 10 elements per wavelength.

## Model limits

- The theory is linear and only valid for small slopes. The CLI prints the peak `|grad u|` and warns above 0.1.
- The pressure wave is a prescribed one-sided load. Air loading and sound radiated by the membrane are ignored, so real drumhead frequencies will differ.
- A true point strike would excite modes the mesh can't resolve, so `strike` uses a finite contact radius and warns if it's too small for the mesh.
- The rim is a polygon inscribed in the circle, which raises frequencies very slightly.
