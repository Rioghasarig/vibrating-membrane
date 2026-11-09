"""Run a membrane scenario from the command line.

Examples:
    python simulate.py strike
    python simulate.py tone --angle 45 --freq 350
    python simulate.py mode --mn 2 1 --view 3d
    python simulate.py blast --save blast.gif --save-summary blast.png --no-show
"""

import argparse
import time

from membrane import simulate
from membrane.scenarios import SCENARIOS


def parse_args():
    parser = argparse.ArgumentParser(
        description="Simulate a vibrating circular membrane with finite elements.",
        epilog="Scenarios: "
        + "; ".join(f"{name}: {fn.__doc__}" for name, fn in SCENARIOS.items()),
    )
    parser.add_argument("scenario", choices=sorted(SCENARIOS))
    parser.add_argument("--rings", type=int, default=40, help="mesh rings (default 40)")
    parser.add_argument("--duration", type=float, help="simulated time [s]")
    parser.add_argument("--dt", type=float, help="time step [s]")
    parser.add_argument("--frames", type=int, help="number of animation frames")
    parser.add_argument("--damping", type=float, help="damping per area [N s/m^3]")
    parser.add_argument("--at", type=float, nargs=2, metavar=("X", "Y"), help="strike point [m]")
    parser.add_argument("--mn", type=int, nargs=2, metavar=("M", "N"), help="Bessel mode indices")
    parser.add_argument("--freq", type=float, help="tone frequency [Hz]")
    parser.add_argument("--pressure", type=float, help="pressure amplitude [Pa]")
    parser.add_argument("--angle", type=float, help="incidence angle from normal [degrees]")
    parser.add_argument("--azimuth", type=float, help="in-plane travel direction [degrees]")
    parser.add_argument("--view", choices=["2d", "3d"], default="2d")
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--save", help="write the animation to this .gif or .mp4")
    parser.add_argument("--save-summary", help="write the summary plots to this image")
    parser.add_argument("--no-show", action="store_true", help="don't open any windows")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.no_show:
        import matplotlib

        matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    from membrane.viz import animate, plot_summary

    setup = SCENARIOS[args.scenario](args)
    m = setup.membrane
    print(
        f"{setup.title}: {m.N} interior nodes, {len(m.triangles)} triangles, "
        f"wave speed {m.wave_speed:.1f} m/s, f01 {m.f01:.1f} Hz"
    )

    start = time.perf_counter()
    result = simulate(
        m,
        q0=setup.q0,
        v0=setup.v0,
        drives=setup.drives,
        duration=setup.duration,
        dt=args.dt,
        n_frames=args.frames or setup.n_frames,
        probes=setup.probes,
    )
    elapsed = time.perf_counter() - start
    print(f"{len(result.step_t) - 1} steps of {result.dt:.3g} s in {elapsed:.2f} s")

    slope = result.max_slope()
    print(f"Peak slope |grad u| = {slope:.3g}")
    if slope > 0.1:
        print("Warning: slope is not small, so the linear theory is unreliable here.")

    # Keep a reference to the animation, or it is garbage collected before it plays.
    anim = None
    if args.save or not args.no_show:
        anim = animate(result, view=args.view, save=args.save, fps=args.fps, title=setup.title)
        if args.save:
            print(f"Wrote {args.save}")
    if args.save_summary or not args.no_show:
        plot_summary(result, title=setup.title, save=args.save_summary)
        if args.save_summary:
            print(f"Wrote {args.save_summary}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
