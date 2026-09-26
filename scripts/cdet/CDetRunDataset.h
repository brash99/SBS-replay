#ifndef CDET_RUN_DATASET_H
#define CDET_RUN_DATASET_H

#include <TChain.h>
#include <TList.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace CDetRunDataset {

struct FileNameInfo {
  int run = -1;
  int firstStream = -1;
  int lastStream = -1;
  int firstSegment = -1;
  int lastSegment = -1;
  int firstEvent = -1;
  int eventCount = -1;
  int part = 0;
  bool isPart = false;
  TString stem;
};

inline bool ParseFileName(const TString &name, FileNameInfo &info) {
  int consumed = 0;
  int part = 0;
  const char *formats[] = {
      "cdet_%d_stream%d_%d_seg%d_%d.root%n",
      "cdet_%d_stream_%d_%d_seg%d_%d.root%n"};
  const char *partFormats[] = {
      "cdet_%d_stream%d_%d_seg%d_%d_%d.root%n",
      "cdet_%d_stream_%d_%d_seg%d_%d_%d.root%n"};
  const char *eventRangeFormats[] = {
      "cdet_%d_stream%d_%d_seg%d_%d_firstevent%d_nevent%d.root%n",
      "cdet_%d_stream_%d_%d_seg%d_%d_firstevent%d_nevent%d.root%n"};
  const char *eventRangePartFormats[] = {
      "cdet_%d_stream%d_%d_seg%d_%d_firstevent%d_nevent%d_%d.root%n",
      "cdet_%d_stream_%d_%d_seg%d_%d_firstevent%d_nevent%d_%d.root%n"};

  for (const char *format : formats) {
    FileNameInfo candidate;
    consumed = 0;
    if (std::sscanf(name.Data(), format, &candidate.run,
                    &candidate.firstStream, &candidate.lastStream,
                    &candidate.firstSegment, &candidate.lastSegment,
                    &consumed) == 5 && consumed == name.Length()) {
      candidate.stem = name(0, name.Length() - 5);
      info = candidate;
      return true;
    }
  }

  for (const char *format : eventRangeFormats) {
    FileNameInfo candidate;
    consumed = 0;
    if (std::sscanf(name.Data(), format, &candidate.run,
                    &candidate.firstStream, &candidate.lastStream,
                    &candidate.firstSegment, &candidate.lastSegment,
                    &candidate.firstEvent, &candidate.eventCount,
                    &consumed) == 7 && consumed == name.Length()) {
      candidate.stem = name(0, name.Length() - 5);
      info = candidate;
      return true;
    }
  }

  for (const char *format : partFormats) {
    FileNameInfo candidate;
    consumed = 0;
    part = 0;
    if (std::sscanf(name.Data(), format, &candidate.run,
                    &candidate.firstStream, &candidate.lastStream,
                    &candidate.firstSegment, &candidate.lastSegment, &part,
                    &consumed) == 6 && consumed == name.Length() && part > 0) {
      candidate.part = part;
      candidate.isPart = true;
      const TString suffix = TString::Format("_%d.root", part);
      candidate.stem = name(0, name.Length() - suffix.Length());
      info = candidate;
      return true;
    }
  }

  for (const char *format : eventRangePartFormats) {
    FileNameInfo candidate;
    consumed = 0;
    part = 0;
    if (std::sscanf(name.Data(), format, &candidate.run,
                    &candidate.firstStream, &candidate.lastStream,
                    &candidate.firstSegment, &candidate.lastSegment,
                    &candidate.firstEvent, &candidate.eventCount, &part,
                    &consumed) == 8 && consumed == name.Length() && part > 0) {
      candidate.part = part;
      candidate.isPart = true;
      const TString suffix = TString::Format("_%d.root", part);
      candidate.stem = name(0, name.Length() - suffix.Length());
      info = candidate;
      return true;
    }
  }
  return false;
}

struct Dataset {
  FileNameInfo info;
  TString baseName;
  std::map<int, TString> parts;
};

inline std::vector<TString> Select(int runNumber, const char *directory,
                                  int segmentMin = -1,
                                  int segmentMax = -1) {
  std::map<std::string, Dataset> datasets;
  TSystemDirectory dir(directory, directory);
  TList *files = dir.GetListOfFiles();
  if (!files)
    return {};

  TIter next(files);
  while (TSystemFile *file = static_cast<TSystemFile *>(next())) {
    if (file->IsDirectory())
      continue;
    const TString name = file->GetName();
    FileNameInfo parsed;
    if (!ParseFileName(name, parsed) || parsed.run != runNumber)
      continue;

    Dataset &dataset = datasets[parsed.stem.Data()];
    if (parsed.isPart) {
      dataset.parts[parsed.part] = name;
    } else {
      dataset.info = parsed;
      dataset.baseName = name;
    }
  }

  bool useRange = segmentMin >= 0 || segmentMax >= 0;
  if (useRange) {
    if (segmentMin < 0)
      segmentMin = segmentMax;
    if (segmentMax < 0)
      segmentMax = segmentMin;
    if (segmentMin > segmentMax)
      std::swap(segmentMin, segmentMax);
  }

  Dataset *selected = nullptr;
  int selectedExcess = 0;
  int selectedSpan = -1;
  for (auto &entry : datasets) {
    Dataset &candidate = entry.second;
    if (candidate.baseName.IsNull())
      continue; // Never select orphaned rollover files.

    const int span = candidate.info.lastSegment - candidate.info.firstSegment;
    if (span < 0)
      continue;

    if (useRange) {
      if (candidate.info.firstSegment > segmentMin ||
          candidate.info.lastSegment < segmentMax)
        continue;
      const int excess = (segmentMin - candidate.info.firstSegment) +
                         (candidate.info.lastSegment - segmentMax);
      if (!selected || excess < selectedExcess ||
          (excess == selectedExcess && span > selectedSpan)) {
        selected = &candidate;
        selectedExcess = excess;
        selectedSpan = span;
      }
    } else if (!selected || span > selectedSpan) {
      selected = &candidate;
      selectedSpan = span;
    }
  }

  if (!selected)
    return {};

  TString prefix(directory);
  if (!prefix.EndsWith("/"))
    prefix += "/";
  std::vector<TString> result;
  result.push_back(prefix + selected->baseName);
  int expectedPart = 1;
  for (const auto &part : selected->parts) {
    if (part.first != expectedPart) {
      std::cerr << "[CDet dataset] Missing rollover part " << expectedPart
                << " for " << selected->baseName << "; stopping before part "
                << part.first << ".\n";
      break;
    }
    result.push_back(prefix + part.second);
    ++expectedPart;
  }
  return result;
}

inline int AddToChain(TChain *chain, int runNumber, const char *directory,
                      int segmentMin = -1, int segmentMax = -1) {
  const std::vector<TString> files =
      Select(runNumber, directory, segmentMin, segmentMax);
  if (files.empty()) {
    std::cerr << "[CDet dataset] No coherent dataset found for run "
              << runNumber << " in " << directory << ".\n";
    return 0;
  }

  std::cout << "[CDet dataset] Selected " << files.size()
            << " file(s) for run " << runNumber << ":\n";
  int added = 0;
  for (const TString &file : files) {
    std::cout << "  " << file << "\n";
    added += chain->Add(file);
  }
  return added;
}

} // namespace CDetRunDataset

#endif
