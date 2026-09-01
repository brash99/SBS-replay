import numpy as np
import glob
import re
import os
import matplotlib.pyplot as plt

# -----------------------------
# User settings
# -----------------------------
calib_dir = "calibrationFiles"

bars_to_plot = [92, 75, 32]      # bars as written in the .dat files
NumBars = 168                    # bars 1 through 168 in the .dat file

file_pattern = f"{calib_dir}/CDet_calibration.dat.*"

output_run_pdf = "bar_offsets_vs_run.pdf"

output_nhits_pdf = "bar_offsets_vs_nhits.pdf"


# -----------------------------
# Helper functions
# -----------------------------
def get_run_number(filename):
    """
    Extract run number from filename like:
    CDet_calibration.dat.4735
    """
    base = os.path.basename(filename)
    match = re.search(r"CDet_calibration\.dat\.(\d+)$", base)

    if match:
        return int(match.group(1))

    return None


def read_bar_offsets(filename):
    """
    Reads the [BarOffsets] chunk from one calibration file.

    Returns:
        offsets[bar - 1, 0] = offset
        offsets[bar - 1, 1] = nhits
    """
    offsets = np.full((NumBars, 2), np.nan)

    in_bar_offsets = False
    found_bar_offsets_section = False
    n_bars_read = 0

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Skip blank lines and comments
            if not line or line.startswith("#"):
                continue

            # Detect section headers more robustly
            if line.startswith("[") and line.endswith("]"):
                section_name = line.strip()

                if section_name == "[BarOffsets]":
                    in_bar_offsets = True
                    found_bar_offsets_section = True
                else:
                    in_bar_offsets = False

                continue

            # Read lines only inside [BarOffsets]
            if in_bar_offsets:
                parts = line.split()

                # Expected:
                # bar offset nhits
                #
                # But this also handles older two-column files:
                # bar offset
                if len(parts) >= 2:
                    try:
                        bar = int(parts[0])
                        offset = float(parts[1])

                        if len(parts) >= 3:
                            nhits = float(parts[2])
                        else:
                            nhits = 0.0

                    except ValueError:
                        print(f"Skipping malformed line in {filename}: {line}")
                        continue

                    if 1 <= bar <= NumBars:
                        row = bar - 1

                        offsets[row, 0] = offset
                        offsets[row, 1] = nhits
                        n_bars_read += 1
                    else:
                        print(f"Warning: bar {bar} in {filename} is outside 1..{NumBars}")

    if not found_bar_offsets_section:
        print(f"Warning: did not find [BarOffsets] section in {filename}")

    if n_bars_read == 0:
        print(f"Warning: read 0 bars from [BarOffsets] in {filename}")

    return offsets, n_bars_read


# -----------------------------
# Main loop
# -----------------------------
data_vs_run = {bar: [] for bar in bars_to_plot}
data_vs_nhits = {bar: [] for bar in bars_to_plot}

files = sorted(glob.glob(file_pattern))

if len(files) == 0:
    raise RuntimeError(f"No files found matching pattern: {file_pattern}")

for filename in files:
    run = get_run_number(filename)

    if run is None:
        print(f"Skipping file with unexpected name: {filename}")
        continue

    offsets, n_bars_read = read_bar_offsets(filename)

    for bar in bars_to_plot:
        if bar < 1 or bar > NumBars:
            print(f"Warning: requested bar {bar} is outside 1..{NumBars}")
            continue

        row = bar - 1

        offset = offsets[row, 0]
        nhits = offsets[row, 1]

        if np.isnan(offset) or np.isnan(nhits):
            print(f"Warning: bar {bar} missing from [BarOffsets] in file {filename}")
            continue

        data_vs_run[bar].append((run, offset))
        data_vs_nhits[bar].append((nhits, offset))


# -----------------------------
# Plot 1: offset vs run number
# -----------------------------
plt.figure(figsize=(8, 5))

for bar in bars_to_plot:
    bar_data = sorted(data_vs_run[bar], key=lambda x: x[0])

    if len(bar_data) == 0:
        print(f"No run data found for bar {bar}")
        continue

    runs = [x[0] for x in bar_data]
    offsets_bar = [x[1] for x in bar_data]

    plt.plot(runs, offsets_bar, marker="o", linestyle="-", label=f"Bar {bar}")

plt.xlabel("Run Number")
plt.ylabel("Bar Offset [ns]")
plt.title("CDet Bar Offsets vs Run Number")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(output_run_pdf)

print(f"Saved plot to {output_run_pdf}")


# -----------------------------
# Plot 2: offset vs nhits
# -----------------------------
plt.figure(figsize=(8, 5))

for bar in bars_to_plot:
    bar_data = sorted(data_vs_nhits[bar], key=lambda x: x[0])

    if len(bar_data) == 0:
        print(f"No nhits data found for bar {bar}")
        continue

    nhits = [x[0] for x in bar_data]
    offsets_bar = [x[1] for x in bar_data]

    plt.scatter(nhits, offsets_bar, marker="o", label=f"Bar {bar}")

plt.xlabel("Number of Hits")
plt.ylabel("Bar Offset [ns]")
plt.title("CDet Bar Offsets vs Number of Hits")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(output_nhits_pdf)

print(f"Saved plot to {output_nhits_pdf}")

plt.show()