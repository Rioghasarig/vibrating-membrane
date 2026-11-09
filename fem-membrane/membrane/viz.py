"""Animation and summary plots."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle
from matplotlib.tri import Triangulation

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]

# Blue for downward, red for upward, neutral gray at rest.
DISPLACEMENT_CMAP = LinearSegmentedColormap.from_list(
    "membrane", ["#104281", "#2a78d6", "#f0efec", "#e34948", "#8f1f1f"]
)


def _length_unit(peak):
    """Pick a readable unit: returns (multiplier from metres, label)."""
    for scale, label in [(1.0, "m"), (1e3, "mm"), (1e6, "µm"), (1e9, "nm")]:
        if peak * scale >= 1:
            return scale, label
    return 1e12, "pm"


def _style_axes(ax):
    ax.set_facecolor(SURFACE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(BASELINE)
    ax.tick_params(colors=MUTED, labelsize=9, length=0)
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(INK_SECONDARY)
    ax.yaxis.label.set_color(INK_SECONDARY)


def _legend(ax):
    # Sits above the plot, level with the panel title, so it can't cover data.
    ax.legend(
        frameon=False,
        fontsize=9,
        labelcolor=INK_SECONDARY,
        loc="lower right",
        bbox_to_anchor=(1, 1),
        ncols=3,
        borderaxespad=0,
    )


def animate(result, view="2d", save=None, fps=30, title=None):
    """Play the simulation in slow motion, or write it to ``save``.

    ``view`` is "2d" for a top-down colour map or "3d" for a surface. The 3d
    view rebuilds the surface every frame, so at high ring counts it runs
    slowly live and is better saved to a file. Returns the animation, which
    the caller must keep a reference to until it has played.
    """
    if view not in ("2d", "3d"):
        raise ValueError("view must be '2d' or '3d'")

    m = result.membrane
    full = m.to_full(result.frames)
    peak = np.abs(full).max()
    if peak == 0:
        raise ValueError("nothing to animate: displacement is zero in every frame")
    scale, unit = _length_unit(peak)
    full = full * scale

    # Clip the colour scale so one sharp early peak doesn't wash out later waves.
    limit = np.percentile(np.abs(full), 99.9)
    limit = limit if limit > 0 else peak * scale
    norm = Normalize(-limit, limit)
    triang = Triangulation(m.x, m.y, m.triangles)

    fig = plt.figure(figsize=(7, 6), facecolor=SURFACE)
    if view == "2d":
        ax = fig.add_subplot()
        ax.set_facecolor(SURFACE)
        mesh = ax.tripcolor(
            triang,
            facecolors=full[0][m.triangles].mean(axis=1),
            cmap=DISPLACEMENT_CMAP,
            norm=norm,
            shading="flat",
            edgecolors="none",
        )
        ax.add_patch(Circle((0, 0), m.a, fill=False, color=BASELINE, linewidth=1.5))
        ax.set_aspect("equal")
        ax.margins(0.03)
        ax.set_axis_off()

        def draw(i):
            mesh.set_array(full[i][m.triangles].mean(axis=1))
    else:
        ax = fig.add_subplot(projection="3d")
        ax.set_facecolor(SURFACE)
        ax.set_zlim(-peak * scale, peak * scale)
        ax.set_box_aspect((1, 1, 0.45))
        ax.set_xlabel("x (m)", color=INK_SECONDARY)
        ax.set_ylabel("y (m)", color=INK_SECONDARY)
        ax.set_zlabel(f"u ({unit})", color=INK_SECONDARY)
        surface = [None]

        def draw(i):
            if surface[0] is not None:
                surface[0].remove()
            surface[0] = ax.plot_trisurf(
                triang,
                full[i],
                cmap=DISPLACEMENT_CMAP,
                vmin=-limit,
                vmax=limit,
                linewidth=0,
                antialiased=False,
            )

    colorbar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=DISPLACEMENT_CMAP),
        ax=ax,
        shrink=0.75,
        pad=0.08,
    )
    colorbar.set_label(f"Displacement ({unit})", color=INK_SECONDARY)
    colorbar.ax.tick_params(colors=MUTED, length=0)
    colorbar.outline.set_visible(False)

    heading = fig.suptitle("", color=INK, fontsize=12)

    def update(i):
        draw(i)
        label = f"t = {result.t[i] * 1e3:.2f} ms"
        heading.set_text(f"{title}   {label}" if title else label)
        return []

    update(0)
    anim = animation.FuncAnimation(
        fig, update, frames=len(full), interval=1000 / fps, blit=False, repeat=True
    )

    if save is not None:
        path = Path(save)
        if path.suffix.lower() == ".gif":
            anim.save(path, writer=animation.PillowWriter(fps=fps))
        elif animation.writers.is_available("ffmpeg"):
            anim.save(path, writer="ffmpeg", fps=fps, dpi=120)
        else:
            raise RuntimeError(
                f"Saving {path.suffix} needs ffmpeg, which wasn't found. "
                "Install it, or save to a .gif instead."
            )
    return anim


def _spectrum(signal, dt):
    """One-sided magnitude spectrum of a Hann-windowed signal."""
    windowed = (signal - signal.mean()) * np.hanning(len(signal))
    return np.fft.rfftfreq(len(signal), dt), np.abs(np.fft.rfft(windowed))


def plot_summary(result, title=None, save=None):
    """Probe displacement, energy and spectrum, as three stacked panels.

    The spectrum is taken from the first probe. The dotted lines are the
    exact mode frequencies of an ideal undamped membrane, labelled (m, n).
    """
    m = result.membrane
    t_ms = result.step_t * 1e3
    fig, (ax_u, ax_e, ax_f) = plt.subplots(
        3, 1, figsize=(8, 8.5), facecolor=SURFACE, layout="constrained"
    )
    if title:
        fig.suptitle(title, color=INK, fontsize=12)
    for ax in (ax_u, ax_e, ax_f):
        _style_axes(ax)

    peak = np.abs(result.probe_u).max()
    scale, unit = _length_unit(peak if peak > 0 else 1.0)
    for i, (x, y) in enumerate(result.probes):
        ax_u.plot(
            t_ms,
            result.probe_u[:, i] * scale,
            color=SERIES[i % len(SERIES)],
            linewidth=1.5,
            label=f"({x * 100:.1f}, {y * 100:.1f}) cm",
        )
    ax_u.set_xlabel("Time (ms)")
    ax_u.set_ylabel(f"Displacement ({unit})")
    ax_u.set_title("Probe displacement", loc="left", color=INK, fontsize=10)
    if len(result.probes) > 1:
        _legend(ax_u)

    for values, label, color in [
        (result.energy, "Total", SERIES[0]),
        (result.kinetic, "Kinetic", SERIES[1]),
        (result.potential, "Potential", SERIES[2]),
    ]:
        ax_e.plot(t_ms, values, color=color, linewidth=1.5, label=label)
    ax_e.set_xlabel("Time (ms)")
    ax_e.set_ylabel("Energy (J)")
    ax_e.set_title("Membrane energy", loc="left", color=INK, fontsize=10)
    _legend(ax_e)

    freqs, magnitude = _spectrum(result.probe_u[:, 0], result.dt)
    f_max = min(6 * m.f01, freqs[-1])
    shown = (freqs > 0) & (freqs <= f_max)
    floor = max(magnitude[shown].max(), np.finfo(float).tiny) * 1e-4
    ax_f.semilogy(
        freqs[shown], np.maximum(magnitude[shown], floor), color=SERIES[0], linewidth=1.5
    )
    ax_f.set_xlim(0, f_max)
    ax_f.set_ylim(top=magnitude[shown].max() * 40)  # headroom for the mode labels
    ax_f.set_xlabel("Frequency (Hz)")
    ax_f.set_ylabel("Magnitude (arbitrary)")
    ax_f.set_title("Spectrum of first probe", loc="left", color=INK, fontsize=10)

    modes = sorted(
        (m.analytic_frequency(mm, nn), mm, nn) for mm in range(9) for nn in range(1, 4)
    )
    modes = [mode for mode in modes if mode[0] <= f_max]
    for rank, (f, mm, nn) in enumerate(modes):
        ax_f.axvline(f, color=MUTED, linewidth=0.6, linestyle=":", zorder=0)
        if rank < 8:  # higher modes are too close together to label
            ax_f.text(
                f,
                0.97,
                f"({mm},{nn})",
                transform=ax_f.get_xaxis_transform(),
                color=INK_SECONDARY,
                fontsize=8,
                rotation=90,
                ha="right",
                va="top",
            )

    if save is not None:
        fig.savefig(save, dpi=150, facecolor=SURFACE)
    return fig
