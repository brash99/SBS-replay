#include <TChain.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CDetRunDataset.h"

namespace {

struct Count {
  Long64_t pulses = 0;
  Long64_t events = 0;
};

void PrintCount(const char *label, const Count &count, Long64_t denominator) {
  const double percent = denominator > 0
                             ? 100.0 * count.pulses / denominator
                             : 0.0;
  std::cout << std::left << std::setw(43) << label << std::right
            << std::setw(12) << count.pulses << std::setw(11)
            << std::fixed << std::setprecision(2) << percent << "%"
            << std::setw(14) << count.events << "\n";
}

} // namespace

void Compare_CDet_PulseCandidate_CutFlow(
    int runNumber = 5710, const char *inputDirectory = nullptr,
    const char *outputFile = "CDet_run5710_pulse_candidate_cutflow.tsv") {
  TString directory(inputDirectory ? inputDirectory : "");
  if (directory.IsNull()) {
    const char *outDir = gSystem->Getenv("OUT_DIR");
    if (outDir)
      directory = outDir;
  }
  if (directory.IsNull()) {
    std::cerr << "[CDet cut flow] An input directory or OUT_DIR is required.\n";
    return;
  }

  TChain chain("T");
  if (CDetRunDataset::AddToChain(&chain, runNumber, directory.Data()) <= 0)
    return;

  TTreeReader reader(&chain);
  TTreeReaderArray<Double_t> pulsePMT(reader, "earm.cdet.pulse.pmtnum");
  TTreeReaderArray<Double_t> valid(reader, "earm.cdet.pulse.calib_valid");
  TTreeReaderArray<Double_t> broad(
      reader, "earm.cdet.pulse.broad_quality_pass");
  TTreeReaderArray<Double_t> ecal(reader, "earm.cdet.pulse.ecal_eligible");
  TTreeReaderArray<Double_t> spatial(reader, "earm.cdet.pulse.spatial_pass");
  TTreeReaderArray<Double_t> legacyPMT(reader, "earm.cdet.hit.pmtnum");
  TTreeReaderArray<Double_t> legacyLE(reader, "earm.cdet.hit.tdc_le");
  TTreeReaderArray<Double_t> legacyToT(reader, "earm.cdet.hit.tdc_tot");
  TTreeReaderArray<Double_t> legacyX(reader, "earm.cdet.hit.xhit");
  TTreeReaderArray<Double_t> legacyY(reader, "earm.cdet.hit.yhit");
  TTreeReaderArray<Double_t> legacyZ(reader, "earm.cdet.hit.zhit");
  TTreeReaderArray<Double_t> legacyMultiplicity(reader, "earm.cdet.tdc_mult");
  TTreeReaderValue<Double_t> ecalX(reader, "earm.ecal.x");
  TTreeReaderValue<Double_t> ecalY(reader, "earm.ecal.y");
  TTreeReaderValue<Double_t> ecalTime(reader, "earm.ecal.adctime");

  Count complete, calibrated, broadCumulative, ecalCumulative;
  Count accepted, acceptedTwoLayer, acceptedUniquePMT;
  Count acceptedUniquePMTTwoLayer, acceptedRepresentedByLegacyHit;
  Count acceptedTwoLayerRepresentedByLegacyHit;
  Count legacyHits, legacyHitsRepresentedByAcceptedPulse;
  Count legacyExact;
  Long64_t differenceNoLegacyHit = 0;
  Long64_t differenceLegacyTiming = 0;
  Long64_t differenceLegacyToT = 0;
  Long64_t differenceLegacyMultiplicity = 0;
  Long64_t differenceLegacyX = 0;
  Long64_t differenceLegacyY = 0;
  Long64_t differenceLegacyEventLayers = 0;
  Long64_t newPulsesPassingLegacyDecision = 0;
  Long64_t commonEventPMTs = 0;
  Long64_t newOnlyEventPMTs = 0;
  Long64_t legacyOnlyEventPMTs = 0;
  Long64_t legacyOnlyNoCompletePulse = 0;
  Long64_t legacyOnlyCalibration = 0;
  Long64_t legacyOnlyBroad = 0;
  Long64_t legacyOnlyECal = 0;
  Long64_t legacyOnlySpatial = 0;
  Long64_t legacyOnlyNewEventLayers = 0;
  Long64_t totalEvents = 0;
  Long64_t malformedEvents = 0;

  while (reader.Next()) {
    ++totalEvents;
    const size_t n = pulsePMT.GetSize();
    if (valid.GetSize() != n || broad.GetSize() != n ||
        ecal.GetSize() != n || spatial.GetSize() != n) {
      ++malformedEvents;
      continue;
    }

    complete.pulses += n;
    if (n > 0)
      ++complete.events;
    legacyHits.pulses += legacyPMT.GetSize();
    if (legacyPMT.GetSize() > 0)
      ++legacyHits.events;

    std::unordered_set<int> legacyIDs;
    std::unordered_map<int, size_t> legacyIndex;
    const size_t nLegacy = std::min(
        {legacyPMT.GetSize(), legacyLE.GetSize(), legacyToT.GetSize(),
         legacyX.GetSize(), legacyY.GetSize(), legacyZ.GetSize(),
         legacyMultiplicity.GetSize()});
    std::vector<bool> legacyBasic(nLegacy, false);
    int legacyLayer1 = 0;
    int legacyLayer2 = 0;
    for (size_t i = 0; i < nLegacy; ++i) {
      legacyIDs.insert(static_cast<int>(legacyPMT[i]));
      legacyIndex[static_cast<int>(legacyPMT[i])] = i;
      const int id = static_cast<int>(legacyPMT[i]);
      const int layer = id < 1344 ? 0 : 1;
      const double correctedX = legacyX[i] * 1.08 - 0.03;
      const bool pass = *ecalY > -1.2 && *ecalY < 1.2 &&
          *ecalX > -1.5 && *ecalX < 1.5 && *ecalX != 0.0 &&
          *ecalY != 0.0 && legacyLE[i] * 0.01 >= 0.02 &&
          legacyLE[i] * 0.01 <= 60.0 && legacyToT[i] * 0.01 >= 4.0 &&
          legacyToT[i] * 0.01 <= 30.0 && legacyMultiplicity[i] < 100.0 &&
          std::fabs(correctedX - (*ecalX) * legacyZ[i] / 6.144) <= 0.08 &&
          std::fabs(legacyY[i] - (*ecalY) * legacyZ[i] / 6.144 - 0.10) <=
              0.36 &&
          *ecalTime > 10.0 && *ecalTime < 35.0;
      legacyBasic[i] = pass;
      if (pass) {
        if (layer == 0)
          ++legacyLayer1;
        else
          ++legacyLayer2;
      }
    }
    const bool legacyTwoLayerEvent = legacyLayer1 >= 1 && legacyLayer1 <= 100 &&
                                     legacyLayer2 >= 1 && legacyLayer2 <= 100;
    if (legacyTwoLayerEvent) {
      legacyExact.events++;
      for (bool pass : legacyBasic)
        legacyExact.pulses += pass;
    }

    Long64_t nValid = 0, nBroad = 0, nECal = 0, nAccepted = 0;
    Long64_t nAcceptedRepresented = 0;
    bool hasLayer1 = false, hasLayer2 = false;
    std::unordered_set<int> acceptedIDs;
    std::unordered_set<int> completeIDs;
    std::unordered_set<int> validIDs;
    std::unordered_set<int> broadIDs;
    std::unordered_set<int> ecalIDs;
    std::unordered_set<int> spatialIDs;
    for (size_t i = 0; i < n; ++i) {
      const int id = static_cast<int>(pulsePMT[i]);
      completeIDs.insert(id);
      const bool passValid = valid[i] > 0.5;
      const bool passBroad = passValid && broad[i] > 0.5;
      const bool passECal = passBroad && ecal[i] > 0.5;
      const bool passAll = passECal && spatial[i] > 0.5;
      if (passValid)
        validIDs.insert(id);
      if (passBroad)
        broadIDs.insert(id);
      if (passECal)
        ecalIDs.insert(id);
      if (passAll)
        spatialIDs.insert(id);
      nValid += passValid;
      nBroad += passBroad;
      nECal += passECal;
      if (!passAll)
        continue;

      ++nAccepted;
      acceptedIDs.insert(id);
      hasLayer1 |= id < 1344;
      hasLayer2 |= id >= 1344;
      if (legacyIDs.count(id)) {
        ++acceptedRepresentedByLegacyHit.pulses;
        ++nAcceptedRepresented;
      }
    }

    calibrated.pulses += nValid;
    broadCumulative.pulses += nBroad;
    ecalCumulative.pulses += nECal;
    accepted.pulses += nAccepted;
    acceptedUniquePMT.pulses += acceptedIDs.size();
    if (nValid > 0)
      ++calibrated.events;
    if (nBroad > 0)
      ++broadCumulative.events;
    if (nECal > 0)
      ++ecalCumulative.events;
    if (nAccepted > 0) {
      ++accepted.events;
      ++acceptedUniquePMT.events;
    }
    if (hasLayer1 && hasLayer2) {
      acceptedTwoLayer.pulses += nAccepted;
      acceptedUniquePMTTwoLayer.pulses += acceptedIDs.size();
      acceptedTwoLayerRepresentedByLegacyHit.pulses += nAcceptedRepresented;
      ++acceptedTwoLayer.events;
      ++acceptedUniquePMTTwoLayer.events;
      if (nAcceptedRepresented > 0)
        ++acceptedTwoLayerRepresentedByLegacyHit.events;

      for (size_t i = 0; i < n; ++i) {
        if (!(valid[i] > 0.5 && broad[i] > 0.5 && ecal[i] > 0.5 &&
              spatial[i] > 0.5))
          continue;
        const int id = static_cast<int>(pulsePMT[i]);
        const auto found = legacyIndex.find(id);
        if (found == legacyIndex.end()) {
          ++differenceNoLegacyHit;
          continue;
        }
        const size_t j = found->second;
        const double correctedX = legacyX[j] * 1.08 - 0.03;
        if (!(legacyLE[j] * 0.01 >= 0.02 && legacyLE[j] * 0.01 <= 60.0))
          ++differenceLegacyTiming;
        else if (!(legacyToT[j] * 0.01 >= 4.0 &&
                   legacyToT[j] * 0.01 <= 30.0))
          ++differenceLegacyToT;
        else if (!(legacyMultiplicity[j] < 100.0))
          ++differenceLegacyMultiplicity;
        else if (!(std::fabs(correctedX - (*ecalX) * legacyZ[j] / 6.144) <=
                   0.08))
          ++differenceLegacyX;
        else if (!(std::fabs(legacyY[j] - (*ecalY) * legacyZ[j] / 6.144 -
                            0.10) <= 0.36))
          ++differenceLegacyY;
        else if (!legacyTwoLayerEvent)
          ++differenceLegacyEventLayers;
        else
          ++newPulsesPassingLegacyDecision;
      }
    }

    std::unordered_set<int> legacyAcceptedIDs;
    if (legacyTwoLayerEvent) {
      for (size_t i = 0; i < nLegacy; ++i) {
        if (legacyBasic[i])
          legacyAcceptedIDs.insert(static_cast<int>(legacyPMT[i]));
      }
    }
    if (hasLayer1 && hasLayer2) {
      for (int id : acceptedIDs) {
        if (legacyAcceptedIDs.count(id))
          ++commonEventPMTs;
        else
          ++newOnlyEventPMTs;
      }
    }
    for (int id : legacyAcceptedIDs) {
      if (hasLayer1 && hasLayer2 && acceptedIDs.count(id))
        continue;
      ++legacyOnlyEventPMTs;
      if (!completeIDs.count(id))
        ++legacyOnlyNoCompletePulse;
      else if (!validIDs.count(id))
        ++legacyOnlyCalibration;
      else if (!broadIDs.count(id))
        ++legacyOnlyBroad;
      else if (!ecalIDs.count(id))
        ++legacyOnlyECal;
      else if (!spatialIDs.count(id))
        ++legacyOnlySpatial;
      else
        ++legacyOnlyNewEventLayers;
    }
    bool anyAcceptedRepresented = false;
    for (int id : acceptedIDs)
      anyAcceptedRepresented |= legacyIDs.count(id) > 0;
    if (anyAcceptedRepresented)
      ++acceptedRepresentedByLegacyHit.events;

    Long64_t representedLegacyThisEvent = 0;
    for (int id : legacyIDs)
      representedLegacyThisEvent += acceptedIDs.count(id) > 0;
    legacyHitsRepresentedByAcceptedPulse.pulses += representedLegacyThisEvent;
    if (representedLegacyThisEvent > 0)
      ++legacyHitsRepresentedByAcceptedPulse.events;
  }

  const Long64_t denominator = complete.pulses;
  std::cout << "\n[CDet pulse-candidate cut flow]\n"
            << "Dataset events: " << totalEvents << "\n"
            << std::left << std::setw(43) << "Population" << std::right
            << std::setw(12) << "entries" << std::setw(12) << "% complete"
            << std::setw(14) << "events" << "\n";
  PrintCount("Complete positive-ToT pulses", complete, denominator);
  PrintCount("+ valid timing calibration", calibrated, denominator);
  PrintCount("+ broad LE/ToT quality", broadCumulative, denominator);
  PrintCount("+ eligible ECal cluster", ecalCumulative, denominator);
  PrintCount("+ projected spatial match (new plot)", accepted, denominator);
  PrintCount("New accepted pulses in two-layer events", acceptedTwoLayer,
             denominator);
  PrintCount("Unique accepted PMTs per event", acceptedUniquePMT, denominator);
  PrintCount("Unique accepted PMTs in two-layer events",
             acceptedUniquePMTTwoLayer, denominator);
  PrintCount("Two-layer accepted pulses with legacy PMT",
             acceptedTwoLayerRepresentedByLegacyHit, denominator);
  PrintCount("Accepted pulses with a legacy hit PMT",
             acceptedRepresentedByLegacyHit, denominator);
  PrintCount("All legacy earm.cdet.hit records", legacyHits, denominator);
  PrintCount("Legacy PMTs with an accepted pulse",
             legacyHitsRepresentedByAcceptedPulse, denominator);
  PrintCount("Exact legacy hAllGoodLe decision", legacyExact, denominator);
  std::cout << "\n[647242 new two-layer pulses classified by legacy decision]\n"
            << "No legacy hit on the same PMT: " << differenceNoLegacyHit << "\n"
            << "Legacy selected hit fails LE window: " << differenceLegacyTiming << "\n"
            << "Legacy selected hit fails ToT window: " << differenceLegacyToT << "\n"
            << "Legacy selected hit fails multiplicity: " << differenceLegacyMultiplicity << "\n"
            << "Legacy selected hit fails x projection: " << differenceLegacyX << "\n"
            << "Legacy selected hit fails y projection: " << differenceLegacyY << "\n"
            << "Legacy event fails two-layer occupancy: " << differenceLegacyEventLayers << "\n"
            << "New pulses passing the legacy decision: " << newPulsesPassingLegacyDecision << "\n";
  std::cout << "\n[Event-PMT population comparison]\n"
            << "Common to new and legacy: " << commonEventPMTs << "\n"
            << "New only: " << newOnlyEventPMTs << "\n"
            << "Legacy only: " << legacyOnlyEventPMTs << "\n"
            << "  no complete pulse: " << legacyOnlyNoCompletePulse << "\n"
            << "  calibration invalid: " << legacyOnlyCalibration << "\n"
            << "  fails broad quality: " << legacyOnlyBroad << "\n"
            << "  fails ECal eligibility: " << legacyOnlyECal << "\n"
            << "  fails new spatial selection: " << legacyOnlySpatial << "\n"
            << "  new accepted set lacks both layers: " << legacyOnlyNewEventLayers << "\n";
  std::cout << "Malformed pulse-array events: " << malformedEvents << "\n";

  std::ofstream out(outputFile);
  out << "population\tentries\tevents\n";
  const auto row = [&out](const char *name, const Count &count) {
    out << name << '\t' << count.pulses << '\t' << count.events << '\n';
  };
  row("complete_positive_tot_pulses", complete);
  row("calibration_valid", calibrated);
  row("broad_quality_cumulative", broadCumulative);
  row("ecal_eligible_cumulative", ecalCumulative);
  row("new_good_pulse_candidates", accepted);
  row("new_candidates_in_two_layer_events", acceptedTwoLayer);
  row("unique_accepted_pmts", acceptedUniquePMT);
  row("unique_accepted_pmts_in_two_layer_events", acceptedUniquePMTTwoLayer);
  row("two_layer_accepted_pulses_with_legacy_pmt",
      acceptedTwoLayerRepresentedByLegacyHit);
  row("accepted_pulses_with_legacy_hit_pmt", acceptedRepresentedByLegacyHit);
  row("legacy_hit_records", legacyHits);
  row("legacy_pmts_with_accepted_pulse", legacyHitsRepresentedByAcceptedPulse);
  row("exact_legacy_hallgoodle_decision", legacyExact);
  out << "difference_no_legacy_hit\t" << differenceNoLegacyHit << "\t0\n"
      << "difference_legacy_le\t" << differenceLegacyTiming << "\t0\n"
      << "difference_legacy_tot\t" << differenceLegacyToT << "\t0\n"
      << "difference_legacy_multiplicity\t" << differenceLegacyMultiplicity << "\t0\n"
      << "difference_legacy_x\t" << differenceLegacyX << "\t0\n"
      << "difference_legacy_y\t" << differenceLegacyY << "\t0\n"
      << "difference_legacy_event_layers\t" << differenceLegacyEventLayers << "\t0\n"
      << "new_pulses_passing_legacy_decision\t" << newPulsesPassingLegacyDecision << "\t0\n";
  out << "common_event_pmts\t" << commonEventPMTs << "\t0\n"
      << "new_only_event_pmts\t" << newOnlyEventPMTs << "\t0\n"
      << "legacy_only_event_pmts\t" << legacyOnlyEventPMTs << "\t0\n"
      << "legacy_only_no_complete_pulse\t" << legacyOnlyNoCompletePulse << "\t0\n"
      << "legacy_only_calibration\t" << legacyOnlyCalibration << "\t0\n"
      << "legacy_only_broad\t" << legacyOnlyBroad << "\t0\n"
      << "legacy_only_ecal\t" << legacyOnlyECal << "\t0\n"
      << "legacy_only_spatial\t" << legacyOnlySpatial << "\t0\n"
      << "legacy_only_new_event_layers\t" << legacyOnlyNewEventLayers << "\t0\n";
  out.close();
  std::cout << "Cut-flow table: " << outputFile << "\n";
}
