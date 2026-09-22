#!/usr/bin/env python3
"""Plot fitted cross-target shift_ns values and the four approved run groups."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cdet_matplotlib_style  # noqa: F401


HERE = Path(__file__).resolve().parent
INPUT = HERE / "CDet_cross_target_shift_calibrations.tsv"
OUTPUT_STEM = HERE / "CDet_cross_target_shift_vs_run"

GROUPS = [
    {
        "name": "Window 1",
        "ecal": r"$-15 < t_{\mathrm{ECal}} < 30$ ns",
        "runs": [3573, 3575],
        "color": "#7b3294",
    },
    {
        "name": "Window 2",
        "ecal": r"$0 < t_{\mathrm{ECal}} < 35$ ns",
        "runs": [
            3602, 3603, 3604, 3648, 3682, 3685, 3686, 3694, 3757,
            3758, 3786, 3788, 3844, 3845, 3846, 4003, 4006,
        ],
        "color": "#008837",
    },
    {
        "name": "Window 3",
        "ecal": r"$10 < t_{\mathrm{ECal}} < 35$ ns",
        "runs": [
            4345, 4400, 4722, 4735, 4926, 4929, 4985, 5066, 5068,
            5290, 5292, 5294, 5423,
        ],
        "color": "#c51b7d",
    },
    {
        "name": "Window 4",
        "ecal": r"$10 < t_{\mathrm{ECal}} < 35$ ns",
        "runs": [5722, 5723, 5724, 5726, 5794, 5795],
        "color": "#0571b0",
    },
]


def load_results() -> dict[int, float]:
    with INPUT.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    failed = [row["run"] for row in rows if row["status"] != "PASS"]
    if failed:
        raise RuntimeError(f"Refusing to plot failed runs: {', '.join(failed)}")

    values = {int(row["run"]): float(row["shift_ns"]) for row in rows}
    expected = {run for group in GROUPS for run in group["runs"]}
    if set(values) != expected:
        missing = sorted(expected - set(values))
        extra = sorted(set(values) - expected)
        raise RuntimeError(f"Run inventory mismatch; missing={missing}, extra={extra}")
    return values


def main() -> None:
    values = load_results()
    all_runs = sorted(values)

    # Boundaries halfway across the gaps between approved run groups make the
    # four acquisition periods visible without implying data in those gaps.
    boundaries = [all_runs[0] - 35]
    for left, right in zip(GROUPS, GROUPS[1:]):
        boundaries.append((left["runs"][-1] + right["runs"][0]) / 2.0)
    boundaries.append(all_runs[-1] + 35)

    fig, ax = plt.subplots(figsize=(12.5, 6.8), constrained_layout=True)

    label_positions = [0.075, 0.27, 0.59, 0.90]
    for index, group in enumerate(GROUPS):
        runs = group["runs"]
        shifts = [values[run] for run in runs]
        mean = sum(shifts) / len(shifts)
        color = group["color"]

        ax.axvspan(
            boundaries[index],
            boundaries[index + 1],
            color=color,
            alpha=0.075,
            linewidth=0,
            zorder=0,
        )
        ax.plot(runs, shifts, color=color, linewidth=1.2, alpha=0.65, zorder=2)
        ax.scatter(
            runs,
            shifts,
            s=45,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            zorder=3,
            label=f'{group["name"]}: {group["ecal"]}',
        )
        ax.hlines(
            mean,
            min(runs),
            max(runs),
            color=color,
            linewidth=2.0,
            linestyle="--",
            zorder=1,
        )
        ax.text(
            label_positions[index],
            0.78,
            f'{group["name"]}\n{group["ecal"]}\n'
            rf'$\langle shift\rangle={mean:.3f}$ ns',
            transform=ax.transAxes,
            color=color,
            fontsize=11.5,
            fontweight="semibold",
            ha="center",
            va="top",
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "white",
                "edgecolor": color,
                "alpha": 0.92,
            },
        )

    for boundary in boundaries[1:-1]:
        ax.axvline(boundary, color="#666666", linewidth=0.8, linestyle=":", zorder=1)

    ax.axhline(0.0, color="black", linewidth=0.8, alpha=0.55)
    ax.set_xlim(boundaries[0], boundaries[-1])
    ax.set_ylim(-13.0, 2.6)
    ax.set_xlabel("Run number", fontsize=14)
    ax.set_ylabel(r"Run timing offset, $shift_{ns}$ (ns)", fontsize=14)
    ax.set_title(
        "CDet cross-target run timing offsets",
        fontsize=17,
        fontweight="semibold",
    )
    ax.text(
        0.01,
        0.02,
        "Dashed segments show the mean within each approved run group.",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#444444",
    )
    ax.grid(axis="y", color="#cccccc", linewidth=0.7, alpha=0.65)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12)
    ax.legend(
        loc="center right",
        fontsize=11,
        frameon=True,
        framealpha=0.95,
        title="Approved ECal ADC-time cuts",
        title_fontsize=11,
    )

    metadata = {
        "Title": "CDet cross-target shift_ns versus run number",
        "Subject": "Four approved cross-target calibration run groups",
        "Creator": "plot_cross_target_shift_vs_run.py",
    }
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=200, metadata=metadata)
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"), metadata=metadata)


if __name__ == "__main__":
    main()
