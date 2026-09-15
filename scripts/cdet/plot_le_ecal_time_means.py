#!/usr/bin/env python3
"""Plot run-average CDet LE and ECal ADC times from the calibration CSV."""

import argparse
import csv


def read_column(filename, value_column):
    values_by_run = {}
    with open(filename, newline="", encoding="utf-8") as csv_file:
        for row in csv.DictReader(csv_file):
            try:
                run = int(row["run_number"])
                value = float(row[value_column])
            except (KeyError, TypeError, ValueError):
                continue
            if value == -999.0:
                continue
            values_by_run[run] = value
    return values_by_run


def main():
    parser = argparse.ArgumentParser(
        description="Plot average CDet LE and ECal ADC time versus run number."
    )
    parser.add_argument(
        "le_csv",
        nargs="?",
        default="cdet_le_means.csv",
        help="CDet LE CSV (default: cdet_le_means.csv)",
    )
    parser.add_argument(
        "ecal_csv",
        nargs="?",
        default="ecal_adctime_means.csv",
        help="ECal ADC-time CSV (default: ecal_adctime_means.csv)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="le_ecal_time_means_vs_run.png",
        help="output image (default: le_ecal_time_means_vs_run.png)",
    )
    args = parser.parse_args()

    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as error:
        raise SystemExit(
            "matplotlib is required; install it with: python3 -m pip install matplotlib"
        ) from error

    le_by_run = read_column(args.le_csv, "mean_le_ns")
    ecal_by_run = read_column(args.ecal_csv, "mean_ecal_adctime_ns")
    runs = sorted(set(le_by_run) & set(ecal_by_run))
    if not runs:
        raise SystemExit("No runs with valid entries in both input CSV files")
    mean_le = [le_by_run[run] for run in runs]
    mean_ecal = [ecal_by_run[run] for run in runs]

    figure, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(runs, mean_le, "o-", markersize=4)
    axes[0].set_ylabel(r"$\langle\mathrm{CDet\ LE}\rangle$ (ns)")
    axes[0].grid(alpha=0.3)

    axes[1].plot(runs, mean_ecal, "o-", markersize=4, color="tab:orange")
    axes[1].set_xlabel("Run number")
    axes[1].set_ylabel(r"$\langle\mathrm{ECal\ ADC\ time}\rangle$ (ns)")
    axes[1].grid(alpha=0.3)

    figure.tight_layout()
    figure.savefig(args.output, dpi=160)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
