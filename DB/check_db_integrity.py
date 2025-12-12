#!/usr/bin/env python3
import sys
import argparse
import re
from typing import List, Dict

SECTIONS = {
    "xpos": "earm.cdet.xpos",
    "ypos": "earm.cdet.ypos",
    "zpos": "earm.cdet.zpos",
}

FLOAT_RE = re.compile(r"""
    (?P<num>
        [+-]?
        (?:
            (?:\d+\.\d*) | (?:\.\d+) | (?:\d+)
        )
        (?:[eE][+-]?\d+)?   # optional exponent
    )
""", re.VERBOSE)


def extract_floats_from_text(s: str) -> List[float]:
    return [float(m.group("num")) for m in FLOAT_RE.finditer(s)]


def parse_section_blankline_terminated(filename: str, section_name: str) -> List[float]:
    """
    Parse a section whose END is indicated by a blank line.
    The section START is a line containing the section header (exact substring match).

    Behavior:
      - Finds the first occurrence of section_name in a line
      - Parses floats appearing on that same line AFTER the header
      - Parses floats on following lines until a blank line (or EOF)
    """
    values: List[float] = []
    in_section = False

    with open(filename, "r") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            stripped = line.strip()

            if not in_section:
                # start when we see the header
                if section_name in stripped:
                    in_section = True
                    # parse any numbers after the header on the same line
                    idx = stripped.find(section_name)
                    tail = stripped[idx + len(section_name):]
                    values.extend(extract_floats_from_text(tail))
                continue

            # in section: blank line ends section
            if stripped == "":
                break

            # otherwise parse numbers on this line
            values.extend(extract_floats_from_text(stripped))

    return values


def main():
    ap = argparse.ArgumentParser(
        description="Verify integrity of CDet position database (blank-line terminated sections)."
    )
    ap.add_argument("dbfile", help="Path to database file (e.g. db_earm.cdet.dat)")
    ap.add_argument("--show", type=int, default=20,
                    help="Max mismatches to print per category (default: 20)")
    ap.add_argument("--tail", type=int, default=16,
                    help="Number of tail entries expected to be +999.000 (default: 16)")
    ap.add_argument("--x_missing", type=float, default=999.0,
                    help="Missing sentinel in xpos core region (default: +999.0)")
    ap.add_argument("--yz_missing", type=float, default=-999.0,
                    help="Missing sentinel in ypos/zpos core region (default: -999.0)")
    args = ap.parse_args()

    data: Dict[str, List[float]] = {}
    for key, sec in SECTIONS.items():
        vals = parse_section_blankline_terminated(args.dbfile, sec)
        if not vals:
            sys.exit(f"ERROR: Could not find or parse any values for section '{sec}'")
        data[key] = vals

    nx, ny, nz = len(data["xpos"]), len(data["ypos"]), len(data["zpos"])

    print("\n=== Entry count check ===")
    print(f"xpos: {nx}")
    print(f"ypos: {ny}")
    print(f"zpos: {nz}")

    if not (nx == ny == nz):
        sys.exit("FAIL: Sections do not have the same number of numerical entries")
    print("PASS: All sections have equal length")

    if nx < args.tail:
        sys.exit(f"ERROR: Section length {nx} is smaller than tail length {args.tail}")

    core_len = nx - args.tail

    # --- Core alignment check ---
    print(f"\n=== Core alignment check (first {core_len} entries) ===")
    mism_y = []
    mism_z = []

    for i in range(core_len):
        x = data["xpos"][i]
        y = data["ypos"][i]
        z = data["zpos"][i]

        if x == args.x_missing:
            if y != args.yz_missing:
                mism_y.append(i)
            if z != args.yz_missing:
                mism_z.append(i)

    if mism_y:
        print(f"FAIL: xpos({args.x_missing:+.3f}) not aligned with ypos({args.yz_missing:+.3f}) "
              f"at {len(mism_y)} indices")
        print(" First indices:", mism_y[:args.show])
    else:
        print(f"PASS: xpos({args.x_missing:+.3f}) aligns with ypos({args.yz_missing:+.3f})")

    if mism_z:
        print(f"FAIL: xpos({args.x_missing:+.3f}) not aligned with zpos({args.yz_missing:+.3f}) "
              f"at {len(mism_z)} indices")
        print(" First indices:", mism_z[:args.show])
    else:
        print(f"PASS: xpos({args.x_missing:+.3f}) aligns with zpos({args.yz_missing:+.3f})")

    # --- Tail check ---
    print(f"\n=== Tail check (last {args.tail} entries must be +999.000) ===")
    for key in ["xpos", "ypos", "zpos"]:
        tail = data[key][-args.tail:]
        bad = [j for j, v in enumerate(tail) if v != 999.0]
        if bad:
            print(f"FAIL: {key} tail not all +999.000")
            # show actual bad values with their absolute indices
            for j in bad[:args.show]:
                abs_i = core_len + j
                print(f"  index {abs_i}: {tail[j]:+.3f}")
            if len(bad) > args.show:
                print(f"  ... ({len(bad) - args.show} more)")
        else:
            print(f"PASS: {key} tail = {args.tail} × +999.000")

    print("\n=== Integrity check complete ===")


if __name__ == "__main__":
    main()
