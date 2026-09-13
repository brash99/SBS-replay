#!/usr/bin/env python3
"""Create the four-run comparison from the lower-left amalgamated panels."""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
RUNS = (
    (5710, "Cross-target", 0.810203, 1.195534,
     "CDet_run5710_shift_diagnostics_bar30_comparison"),
    (5992, "Cross-target", 0.810203, 5.171528,
     "CDet_run5992_shift_diagnostics_bar30_comparison"),
    (5711, "LH2", -0.012423, 3.422000,
     "CDet_run5711_shift_diagnostics"),
    (6077, "LH2", 0.002520, 2.913449,
     "CDet_run6077_shift_diagnostics_final"),
)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    name = ("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else
            "/System/Library/Fonts/Supplemental/Arial.ttf")
    return ImageFont.truetype(name, size)


def lower_left_panel(path: Path) -> Image.Image:
    source = Image.open(path).convert("RGB")
    width, height = source.size
    # ROOT divides this canvas into equal 2x2 pads. Include a small strip above
    # the midpoint so the complete lower-left panel title is retained.
    panel = source.crop((0, int(height * 0.485), int(width * 0.505), height))
    return panel.resize((1200, 690), Image.Resampling.LANCZOS)


def main() -> None:
    canvas = Image.new("RGB", (2500, 1700), "white")
    draw = ImageDraw.Draw(canvas)
    draw.text(
        (1250, 34),
        "CDet Bar 30: fully selected ECal–CDet timing spectra",
        anchor="ma",
        font=font(46, True),
        fill="black",
    )

    positions = ((35, 135), (1265, 135), (35, 905), (1265, 905))
    for (run, target, p1, shift, directory), (x, y) in zip(RUNS, positions):
        source = (HERE / directory / "bar_ecal_cdet_timing" /
                  f"run_{run}_stage_7" / "bar_030_amalgamated.png")
        panel = lower_left_panel(source)
        canvas.paste(panel, (x, y + 82))
        draw.text(
            (x + 600, y), f"Run {run} — {target}", anchor="ma",
            font=font(35, True), fill="black",
        )
        draw.text(
            (x + 600, y + 43),
            f"p1 = {p1:.6f} ns/ns     shift_ns = {shift:.6f} ns",
            anchor="ma", font=font(27), fill="black",
        )
        draw.rectangle((x, y + 82, x + 1200, y + 772), outline="black", width=2)

    png = HERE / "CDet_bar30_four_run_comparison.png"
    pdf = HERE / "CDet_bar30_four_run_comparison.pdf"
    canvas.save(png, dpi=(180, 180))
    canvas.save(pdf, "PDF", resolution=180.0)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
