#!/usr/bin/env python3
"""Compare representative LH2 timing spectra with their CDet era fallbacks."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

import cdet_matplotlib_style  # noqa: F401


HERE = Path(__file__).resolve().parent
HISTOGRAM = "hCDetSelectedPairMeanCorrectedLE"
REFERENCE_RUN = 5711
REFERENCE_SHIFT_NS = 3.422000
OUTPUT_STEM = HERE / "CDet_lh2_shift_relative_to_era"
OUTPUT_TABLE = OUTPUT_STEM.with_suffix(".tsv")

# The replay shift is the value actually applied when each survey ROOT file was
# produced. All survey runs except the independently calibrated Run 5711 used
# their contemporaneous cross-target era mean.
RUNS = {
    3649: {"era": 2, "era_shift": -9.359456, "replay_shift": -9.359456,
           "calibrated": False},
    4346: {"era": 3, "era_shift": -11.060945, "replay_shift": -11.060945,
           "calibrated": False},
    4724: {"era": 3, "era_shift": -11.060945, "replay_shift": -11.060945,
           "calibrated": False},
    5295: {"era": 3, "era_shift": -11.060945, "replay_shift": -11.060945,
           "calibrated": False},
    5711: {"era": 4, "era_shift": 1.095437, "replay_shift": REFERENCE_SHIFT_NS,
           "calibrated": True},
    5727: {"era": 4, "era_shift": 1.095437, "replay_shift": 1.095437,
           "calibrated": False},
    5886: {"era": 4, "era_shift": 1.095437, "replay_shift": 1.095437,
           "calibrated": False},
    6077: {"era": 4, "era_shift": 1.095437, "replay_shift": 2.913449,
           "calibrated": True},
}

ERA_COLORS = {2: "#008837", 3: "#c51b7d", 4: "#0571b0"}


def read_histogram(run: int) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    path = HERE / f"CDet_run{run}_good_pulse_tdc" / "CDetGoodPulse_AllTDC.root"
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    histogram = root_file.Get(HISTOGRAM)
    if not histogram:
        raise RuntimeError(f"{HISTOGRAM} is missing from {path}")

    bins = histogram.GetNbinsX()
    centers = np.array(
        [histogram.GetBinCenter(index) for index in range(1, bins + 1)]
    )
    counts = np.array(
        [histogram.GetBinContent(index) for index in range(1, bins + 1)]
    )
    probability = np.array([0.5], dtype="d")
    quantile = np.zeros(1, dtype="d")
    histogram.GetQuantiles(1, quantile, probability)
    entries = float(histogram.GetEntries())
    mean = float(histogram.GetMean())
    median = float(quantile[0])
    root_file.Close()
    return centers, counts, entries, mean, median


def main() -> None:
    spectra: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    results: dict[int, dict[str, float]] = {}
    for run, settings in RUNS.items():
        centers, counts, entries, mean, median = read_histogram(run)
        spectra[run] = (centers, counts)
        results[run] = {
            **settings,
            "entries": entries,
            "mean": mean,
            "median": median,
        }

    reference_median = results[REFERENCE_RUN]["median"]
    for run, result in results.items():
        median_correction = reference_median - result["median"]
        inferred_shift = result["replay_shift"] + median_correction
        result["median_correction"] = median_correction
        result["inferred_shift"] = inferred_shift
        result["comparison_shift"] = (
            result["replay_shift"] if result["calibrated"] else inferred_shift
        )
        result["increment_above_era"] = (
            result["comparison_shift"] - result["era_shift"]
        )

    with OUTPUT_TABLE.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "run",
            "era",
            "independently_calibrated",
            "entries",
            "era_shift_ns",
            "replay_shift_ns",
            "pair_time_mean_ns",
            "pair_time_median_ns",
            "median_correction_ns",
            "median_aligned_shift_ns",
            "comparison_shift_ns",
            "increment_above_era_ns",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for run, result in results.items():
            writer.writerow(
                {
                    "run": run,
                    "era": int(result["era"]),
                    "independently_calibrated": int(result["calibrated"]),
                    "entries": int(result["entries"]),
                    "era_shift_ns": f'{result["era_shift"]:.6f}',
                    "replay_shift_ns": f'{result["replay_shift"]:.6f}',
                    "pair_time_mean_ns": f'{result["mean"]:.6f}',
                    "pair_time_median_ns": f'{result["median"]:.6f}',
                    "median_correction_ns": f'{result["median_correction"]:.6f}',
                    "median_aligned_shift_ns": f'{result["inferred_shift"]:.6f}',
                    "comparison_shift_ns": f'{result["comparison_shift"]:.6f}',
                    "increment_above_era_ns":
                        f'{result["increment_above_era"]:.6f}',
                }
            )

    fig = plt.figure(figsize=(14.5, 9.0), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.2, 1.0))
    spectra_axis = fig.add_subplot(grid[0, :])
    shift_axis = fig.add_subplot(grid[1, 0])
    increment_axis = fig.add_subplot(grid[1, 1])

    linestyles = ["-", "--", "-.", ":", "-", "--", ":", "-."]
    for (run, result), linestyle in zip(results.items(), linestyles):
        centers, counts = spectra[run]
        normalized = counts / counts.sum()
        color = "#222222" if run == REFERENCE_RUN else ERA_COLORS[int(result["era"])]
        width = 2.8 if run == REFERENCE_RUN else 1.8
        spectra_axis.step(
            centers,
            normalized,
            where="mid",
            color=color,
            linewidth=width,
            linestyle=linestyle,
            label=(
                f'Run {run}, Era {int(result["era"])}: '
                f'median {result["median"]:.2f} ns'
            ),
        )
    spectra_axis.axvline(30.0, color="#555555", linestyle=":", linewidth=1.5)
    spectra_axis.text(30.25, 0.97, "30 ns reference", rotation=90,
                      transform=spectra_axis.get_xaxis_transform(), va="top",
                      color="#555555")
    spectra_axis.set_xlim(8.0, 48.0)
    spectra_axis.set_xlabel(r"Absolute corrected selected-pair mean time (ns)")
    spectra_axis.set_ylabel("Normalized pairs per 1 ns")
    spectra_axis.set_title("Representative LH2 timing spectra with the current database constants")
    spectra_axis.legend(ncol=2, fontsize=10.2, framealpha=0.95)
    spectra_axis.grid(alpha=0.25)

    runs = np.array(list(results))
    era_shifts = np.array([results[run]["era_shift"] for run in runs])
    comparison_shifts = np.array(
        [results[run]["comparison_shift"] for run in runs]
    )
    increments = np.array([results[run]["increment_above_era"] for run in runs])
    colors = [ERA_COLORS[int(results[run]["era"])] for run in runs]

    for run, baseline, comparison, color in zip(
        runs, era_shifts, comparison_shifts, colors
    ):
        shift_axis.plot([run, run], [baseline, comparison], color=color,
                        linewidth=2.3, alpha=0.8)
    shift_axis.scatter(runs, era_shifts, marker="s", s=60, color=colors,
                       edgecolor="white", linewidth=0.7,
                       label="Established era fallback")
    shift_axis.scatter(runs, comparison_shifts, marker="o", s=75, color=colors,
                       edgecolor="black", linewidth=0.6,
                       label="Independently calibrated or Run-5711-inferred LH2 shift")
    shift_axis.set_xlabel("Run number")
    shift_axis.set_ylabel(r"$shift_{ns}$ (ns)")
    shift_axis.set_title("Era fallback and inferred LH2 shift")
    shift_axis.grid(alpha=0.25)
    shift_axis.legend(fontsize=9.5)

    increment_axis.bar([str(run) for run in runs], increments, color=colors,
                       edgecolor="white", linewidth=0.8)
    increment_axis.axhline(0.0, color="black", linewidth=0.8)
    for index, value in enumerate(increments):
        increment_axis.text(index, value + 0.10, f"{value:.2f}",
                            ha="center", va="bottom", fontsize=10)
    increment_axis.set_ylim(0.0, max(increments) + 0.8)
    increment_axis.set_xlabel("LH2 survey run")
    increment_axis.set_ylabel(r"Additional shift above era fallback (ns)")
    increment_axis.set_title("Quantity of interest: LH2 increment relative to its era")
    increment_axis.grid(axis="y", alpha=0.25)

    fig.suptitle(
        "CDet LH2 timing relative to the four cross-target calibration eras",
        fontsize=17,
        fontweight="semibold",
    )
    metadata = {
        "Title": "CDet LH2 shift relative to era fallback",
        "Subject": "Representative LH2 timing spectra and inferred shift increments",
        "Creator": "plot_cdet_lh2_shift_relative_to_era.py",
    }
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=200, metadata=metadata)
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"), metadata=metadata)


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    main()
