import numpy as np
import glob
import re
import os
import matplotlib.pyplot as plt

# -----------------------------
# User settings
# -----------------------------
calib_dir = "calibrationFiles"

bars_to_plot = [92, 75, 30]   # bars as written in the .dat files
NumBars = 168                 # .dat files are 1-indexed: bars 1 through 168

file_pattern = f"{calib_dir}/CDet_calibration.dat.*"


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

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Skip blank lines and comments
            if not line or line.startswith("#"):
                continue

            # Detect section headers
            if line.startswith("[") and line.endswith("]"):
                in_bar_offsets = (line == "[BarOffsets]")
                continue

            # Read lines only inside [BarOffsets]
            if in_bar_offsets:
                parts = line.split()

                # Supports either:
                # bar offset nhits
                # or older files:
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
                    else:
                        print(f"Warning: bar {bar} in {filename} is outside 1..{NumBars}")

    return offsets


# -----------------------------
# Main loop
# -----------------------------
# For each bar, store tuples of (run, offset, nhits)
data = {bar: [] for bar in bars_to_plot}

files = sorted(glob.glob(file_pattern))

if len(files) == 0:
    raise RuntimeError(f"No files found matching pattern: {file_pattern}")

for filename in files:
    run = get_run_number(filename)

    if run is None:
        print(f"Skipping file with unexpected name: {filename}")
        continue

    offsets = read_bar_offsets(filename)

    for bar in bars_to_plot:
        if bar < 1 or bar > NumBars:
            print(f"Warning: requested bar {bar} is outside 1..{NumBars}")
            continue

        row = bar - 1

        offset = offsets[row, 0]
        nhits = offsets[row, 1]

        # Missing values remain np.nan
        # nhits = 0 is still valid and should be plotted
        if np.isnan(offset) or np.isnan(nhits):
            print(f"Warning: bar {bar} missing from [BarOffsets] in file {filename}")
            continue

        data[bar].append((run, offset, nhits))


# -----------------------------
# Plot one heatmap-style figure per bar
# -----------------------------
for bar in bars_to_plot:
    bar_data = sorted(data[bar], key=lambda x: x[0])

    if len(bar_data) == 0:
        print(f"No data found for bar {bar}")
        continue

    runs = [x[0] for x in bar_data]
    offsets_bar = [x[1] for x in bar_data]
    nhits_bar = [x[2] for x in bar_data]

    plt.figure(figsize=(8, 5))

    sc = plt.scatter(runs, offsets_bar, c=nhits_bar, s=60)

    plt.xlabel("Run Number")
    plt.ylabel("Bar Offset [ns]")
    plt.title(f"Bar {bar}: Offset vs Run Number (color = nhits)")
    plt.grid(True, alpha=0.3)

    cbar = plt.colorbar(sc)
    cbar.set_label("Number of Hits")

    plt.tight_layout()

    output_pdf = f"bar_{bar}_offset_vs_run_nhits.pdf"

    plt.savefig(output_pdf)
    plt.savefig(output_png)

    print(f"Saved plot for bar {bar} to {output_pdf} and {output_png}")

plt.show()