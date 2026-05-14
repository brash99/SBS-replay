import glob
import re
import matplotlib.pyplot as plt

# -----------------------------
# User settings
# -----------------------------
calib_dir = "calibrationFiles"   # directory containing CDet_calibration.dat.runnum files

bars_to_plot = [92, 75, 32]   # change these to whichever 3 bars you want

file_pattern = f"{calib_dir}/CDet_calibration.dat.*"

output_pdf = "bar_offsets_vs_run.pdf"
output_png = "bar_offsets_vs_run.png"


# -----------------------------
# Helper functions
# -----------------------------
def get_run_number(filename):
    """
    Extract run number from filename like:
    CDet_calibration.dat.4735
    """
    match = re.search(r"CDet_calibration\.dat\.(\d+)$", filename)
    if match:
        return int(match.group(1))
    return None


def read_bar_offsets(filename):
    """
    Reads the [BarOffsets] chunk from one calibration file.

    Returns:
        dict: {bar_number: offset}
    """
    offsets = {}
    in_bar_offsets = False

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # skip blank lines and comments
            if not line or line.startswith("#"):
                continue

            # detect section headers
            if line.startswith("[") and line.endswith("]"):
                if line == "[BarOffsets]":
                    in_bar_offsets = True
                else:
                    in_bar_offsets = False
                continue

            # read lines only inside [BarOffsets]
            if in_bar_offsets:
                parts = line.split()

                if len(parts) >= 2:
                    bar = int(parts[0])
                    offset = float(parts[1])
                    offsets[bar] = offset

    return offsets


# -----------------------------
# Main loop
# -----------------------------
data = {bar: [] for bar in bars_to_plot}

files = glob.glob(file_pattern)

if len(files) == 0:
    raise RuntimeError(f"No files found matching pattern: {file_pattern}")

for filename in files:
    run = get_run_number(filename)

    if run is None:
        print(f"Skipping file with unexpected name: {filename}")
        continue

    offsets = read_bar_offsets(filename)

    for bar in bars_to_plot:
        if bar in offsets:
            data[bar].append((run, offsets[bar]))
        else:
            print(f"Warning: bar {bar} not found in file {filename}")


# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(8, 5))

for bar in bars_to_plot:
    # sort by run number
    bar_data = sorted(data[bar], key=lambda x: x[0])

    if len(bar_data) == 0:
        print(f"No data found for bar {bar}")
        continue

    runs = [x[0] for x in bar_data]
    offsets = [x[1] for x in bar_data]

    plt.scatter(runs, offsets, marker="o", linestyle="-", label=f"Bar {bar}")

plt.xlabel("Run Number")
plt.ylabel("Bar Offset [ns]")
plt.title("CDet Bar Offsets vs Run Number")
plt.grid(True, alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig(output_pdf)
plt.savefig(output_png)

plt.show()

print(f"Saved plot to {output_pdf} and {output_png}")