#ifndef CDET_GOOD_PULSE_CONFIG_H
#define CDET_GOOD_PULSE_CONFIG_H

#include <TEnv.h>

#include <iostream>

namespace CDetGoodPulseConfig {

struct Values {
  int runNumber = 5711;
  double binWidthNs = 1.0;
  double leMinNs = 0.0;
  double leMaxNs = 60.0;
  double totMinNs = 0.0;
  double totMaxNs = 40.0;
  bool recoveredOnly = false;
  double pairResidualCenterM = 0.0;
  double pairTimingCenterNs = -26.0;
  double pairResidualScaleM = 0.020;
  double pairTimingScaleNs = 5.0;
  double pairCutRadius = 2.0;
  double singleResidualCenterM = 0.0;
  double singleTimingCenterNs = -26.0;
  double singleResidualScaleM = 0.040;
  double singleTimingScaleNs = 5.0;
  double singleCutRadius = 2.0;
  bool oppositeSideEnabled = false;
  double oppositeDYCenterM = 0.51;
  double oppositeDYToleranceM = 0.08;
  double oppositeProjectedYCenterM = 0.0;
  double oppositeProjectedYMaxM = 0.17;
  double ecalTimeMinNs = -10.0;
  double ecalTimeMaxNs = 10.0;
  double ecalEnergyMinGeV = 3.0;
  double ecalEnergyMaxGeV = 4.5;
  double scanMinHalfWidthNs = 0.5;
  double scanMaxHalfWidthNs = 10.0;
  double scanStepNs = 0.5;
};

inline bool Load(const char *configFile, Values &values,
                 const char *caller = "CDet good-pulse configuration") {
  if (!configFile || !configFile[0]) {
    std::cerr << '[' << caller << "] ERROR: configuration filename is empty.\n";
    return false;
  }

  TEnv env;
  if (env.ReadFile(configFile, kEnvLocal) != 0) {
    std::cerr << '[' << caller << "] ERROR: cannot read configuration file "
              << configFile << ".\n";
    return false;
  }
  if (env.GetValue("config.version", 0) != 1) {
    std::cerr << '[' << caller << "] ERROR: " << configFile
              << " must set config.version = 1.\n";
    return false;
  }

  values.runNumber = env.GetValue("analysis.run_number", values.runNumber);
  values.ecalTimeMinNs =
      env.GetValue("analysis.ecal_time_min", values.ecalTimeMinNs);
  values.ecalTimeMaxNs =
      env.GetValue("analysis.ecal_time_max", values.ecalTimeMaxNs);
  values.ecalEnergyMinGeV =
      env.GetValue("analysis.ecal_energy_min", values.ecalEnergyMinGeV);
  values.ecalEnergyMaxGeV =
      env.GetValue("analysis.ecal_energy_max", values.ecalEnergyMaxGeV);

  values.binWidthNs =
      env.GetValue("good_pulse.bin_width_ns", values.binWidthNs);
  values.leMinNs = env.GetValue("good_pulse.le_min_ns", values.leMinNs);
  values.leMaxNs = env.GetValue("good_pulse.le_max_ns", values.leMaxNs);
  values.totMinNs = env.GetValue("good_pulse.tot_min_ns", values.totMinNs);
  values.totMaxNs = env.GetValue("good_pulse.tot_max_ns", values.totMaxNs);
  values.recoveredOnly =
      env.GetValue("good_pulse.recovered_only", values.recoveredOnly ? 1 : 0) != 0;
  values.pairResidualCenterM = env.GetValue(
      "good_pulse.pair_residual_center_m", values.pairResidualCenterM);
  values.pairTimingCenterNs = env.GetValue(
      "good_pulse.pair_timing_center_ns", values.pairTimingCenterNs);
  values.pairResidualScaleM = env.GetValue(
      "good_pulse.pair_residual_scale_m", values.pairResidualScaleM);
  values.pairTimingScaleNs = env.GetValue(
      "good_pulse.pair_timing_scale_ns", values.pairTimingScaleNs);
  values.pairCutRadius = env.GetValue(
      "good_pulse.pair_cut_radius", values.pairCutRadius);
  values.singleResidualCenterM = env.GetValue(
      "good_pulse.single_residual_center_m", values.singleResidualCenterM);
  values.singleTimingCenterNs = env.GetValue(
      "good_pulse.single_timing_center_ns", values.singleTimingCenterNs);
  values.singleResidualScaleM = env.GetValue(
      "good_pulse.single_residual_scale_m", values.singleResidualScaleM);
  values.singleTimingScaleNs = env.GetValue(
      "good_pulse.single_timing_scale_ns", values.singleTimingScaleNs);
  values.singleCutRadius = env.GetValue(
      "good_pulse.single_cut_radius", values.singleCutRadius);
  values.oppositeSideEnabled = env.GetValue(
      "good_pulse.opposite_side_enable",
      values.oppositeSideEnabled ? 1 : 0) != 0;
  values.oppositeDYCenterM = env.GetValue(
      "good_pulse.opposite_dy_center_m", values.oppositeDYCenterM);
  values.oppositeDYToleranceM = env.GetValue(
      "good_pulse.opposite_dy_tolerance_m", values.oppositeDYToleranceM);
  values.oppositeProjectedYCenterM = env.GetValue(
      "good_pulse.opposite_projected_y_center_m",
      values.oppositeProjectedYCenterM);
  values.oppositeProjectedYMaxM = env.GetValue(
      "good_pulse.opposite_projected_y_max_m",
      values.oppositeProjectedYMaxM);

  values.scanMinHalfWidthNs = env.GetValue(
      "pair_scan.min_half_width_ns", values.scanMinHalfWidthNs);
  values.scanMaxHalfWidthNs = env.GetValue(
      "pair_scan.max_half_width_ns", values.scanMaxHalfWidthNs);
  values.scanStepNs =
      env.GetValue("pair_scan.step_ns", values.scanStepNs);

  const bool valid = values.runNumber > 0 && values.binWidthNs > 0.0 &&
      values.leMaxNs > values.leMinNs && values.totMaxNs > values.totMinNs &&
      values.pairResidualScaleM > 0.0 && values.pairTimingScaleNs > 0.0 &&
      values.pairCutRadius > 0.0 && values.singleResidualScaleM > 0.0 &&
      values.singleTimingScaleNs > 0.0 && values.singleCutRadius > 0.0 &&
      (!values.oppositeSideEnabled ||
       (values.oppositeDYCenterM > 0.0 &&
        values.oppositeDYToleranceM > 0.0 &&
        values.oppositeProjectedYMaxM > 0.0)) &&
      values.ecalTimeMaxNs > values.ecalTimeMinNs &&
      values.ecalEnergyMaxGeV > values.ecalEnergyMinGeV &&
      values.scanMinHalfWidthNs > 0.0 &&
      values.scanMaxHalfWidthNs >= values.scanMinHalfWidthNs &&
      values.scanStepNs > 0.0;
  if (!valid) {
    std::cerr << '[' << caller << "] ERROR: invalid good-pulse or pair-scan "
              << "parameters in " << configFile << ".\n";
    return false;
  }
  return true;
}

} // namespace CDetGoodPulseConfig

#endif
