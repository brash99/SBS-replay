#pragma once

#include <TStyle.h>

// Shared presentation defaults for CDet ROOT plots. ROOT's stock axis text is
// too small on the multi-pad canvases used by these analysis macros.
inline void ApplyCDetPlotStyle()
{
  if (!gStyle) return;

  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");

  // Sizes are fractions of the pad dimensions. These remain legible both in
  // an interactive canvas and after exporting the canvas to PDF or PNG.
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetTitleOffset(1.05, "X");
  gStyle->SetTitleOffset(1.25, "Y");
  gStyle->SetTitleOffset(1.20, "Z");

  // Leave enough room for the larger labels and axis titles.
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadTopMargin(0.09);

  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.040);
  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.035);
}

namespace {
struct CDetPlotStyleInitializer {
  CDetPlotStyleInitializer() { ApplyCDetPlotStyle(); }
};

const CDetPlotStyleInitializer gCDetPlotStyleInitializer;
}

