// Let's compare elastic yields for different GEM HV settings:
//
// Run 2977 GEM HV at 3550 (85%)
// Run 2980 GEM HV at 3575 except L13M23 L13M1 L15M01
//
// Alleged elastic Cut:
// root [4] C->Draw("heep.dxECAL","(((((abs(heep.dxECAL-0.03)<0.06 && abs(heep.dyECAL-0.023)<0.08) && (abs(heep.dt_ADC-115)<20.0 && abs(heep.dpp+0.12)<0.1)) && (sbs.gemFT.track.chi2ndf[0]<30)) && ((sbs.gemFT.track.nhits[0]>4 || sbs.gemFT.track.chi2ndf[0]<15))) && (abs(sbs.tr.x[0] + sbs.tr.th[0]*sbs.z_bcp[30] - sbs.x_bcp[30] + 0.076)<0.21)) && (abs(sbs.tr.y[0] + sbs.tr.ph[0]*sbs.z_bcp[30] - sbs.y_bcp[30])<0.21)");
