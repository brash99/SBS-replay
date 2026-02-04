#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>      // for sscanf
#include <algorithm>   // for std::sort
#include <TMath.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TLatex.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>


using namespace std;

const TString REPLAYED_DIR = TString(gSystem->Getenv("OUT_DIR")) + "/rootfiles";

// const TString ANALYSED_DIR = gSystem->Getenv("ANALYSED_DIR");
//const TString REPLAYED_DIR = "/volatile/halla/sbs/btspaude/cdet/rootfiles";
const TString ANALYSED_DIR = "/work/halla/sbs/btspaude/sbs/Rootfiles/cdetFiles";

// Parse the "segX_Y" part: returns true and fills firstSeg/lastSeg if found.
bool GetSegRange(const TString& fname, int& firstSeg, int& lastSeg) {
  // Find "_seg"
  Ssiz_t pos = fname.Index("_seg");
  if (pos == kNPOS) return false;

  // Tail looks like "9_9.root" or "9_9_1.root"
  TString tail = fname(pos + 4, fname.Length() - (pos + 4));

  // Extract first two ints; ignore any further suffix
  int a = -1, b = -1;
  if (sscanf(tail.Data(), "%d_%d", &a, &b) == 2) {
    firstSeg = a; lastSeg = b;
    return true;
  }
  return false;
}

void AddRunFilesToChain(TChain *chain, const char *dir, int runnum, int segMin = -1, int segMax = -1) {
  TString prefix = dir;
  std::vector<TString> runfiles;

  TSystemDirectory directory(prefix, prefix);
  TList *files = directory.GetListOfFiles();

  if (files) {
    TIter next(files);
    TSystemFile *f;
    while ((f = (TSystemFile*) next())) {
      if (f->IsDirectory()) continue; // skip dirs like "." and ".."

      TString fname = f->GetName();
      if (!fname.BeginsWith(Form("cdet_%d_", runnum))) continue;
      if (!fname.EndsWith(".root")) continue;

      // Range filtering enabled only if segMin/segMax are set
      if (segMin >= 0 || segMax >= 0) {
        if (segMin < 0) segMin = segMax;
        if (segMax < 0) segMax = segMin;
        if (segMin > segMax) std::swap(segMin, segMax);

        int firstSeg = -1, lastSeg = -1;
        if (!GetSegRange(fname, firstSeg, lastSeg)) continue;

        // accept if [firstSeg,lastSeg] overlaps [segMin,segMax]
        if (lastSeg < segMin || firstSeg > segMax) continue;
      }

      runfiles.push_back(prefix + "/" + fname);
    }
  }

  std::sort(runfiles.begin(), runfiles.end());

  std::cout << "Adding " << runfiles.size() << " files for run " << runnum << "...\n";
  for (auto &file : runfiles) {
    std::cout << "  " << file << "\n";
    chain->Add(file);
  }
}

std::vector<int> getLocation(int pixelID) {
  // Check valid range
  if (pixelID < 0 || pixelID > 2687) {
      std::cerr << "Error: pixelID must be in the range 0 to 2687.\n";
      return {};  // return empty vector to signal error
  }

  int layerNum      = pixelID / 1344; //1344 pixels per layer
  pixelID               %= 1344;

  int sideNum       = pixelID / 672;  //672 pixels per side
  pixelID               %= 672;

  int submoduleNum  = pixelID / 224; //224 pixels per side of module
  pixelID               %= 224;

  int pmtNum        = pixelID / 16; //16 pixels per bar
  pixelID               %= 16;
  int pixelNum      = pixelID % 16;

  return {layerNum, sideNum, submoduleNum, pmtNum, pixelNum};
}

static const std::unordered_set<int> kUnusedCDetPixels = {
3, 13, 28, 31, 41, 42, 57, 59, 65, 79, 83, 95, 109, 111, 115, 127,
140, 143, 145, 156, 172, 175, 176, 188, 195, 199, 213, 220, 236, 239, 244, 255,
268, 271, 284, 287, 300, 303, 307, 319, 332, 335, 339, 351, 354, 364, 371, 381,
384, 396, 401, 410, 419, 423, 435, 436, 451, 461, 465, 479, 480, 483, 508, 511,
512, 515, 540, 543, 546, 559, 563, 573, 576, 589, 596, 605, 609, 610, 627, 638,
643, 655, 656, 665, 674, 675, 696, 703, 707, 709, 725, 729, 738, 748, 752, 766,
777, 780, 784, 791, 800, 812, 818, 828, 844, 847, 850, 860, 867, 868, 884, 885,
900, 904, 912, 927, 940, 943, 945, 947, 967, 971, 986, 991, 1005, 1007, 1011, 1023,
1027, 1028, 1043, 1050, 1066, 1068, 1072, 1075, 1088, 1102, 1106, 1119, 1121, 1135, 1148, 1151,
1162, 1166, 1178, 1182, 1184, 1186, 1203, 1215, 1228, 1231, 1235, 1247, 1249, 1252, 1267, 1274,
1282, 1285, 1299, 1310, 1317, 1321, 1340, 1341, 1349, 1359, 1372, 1375, 1376, 1391, 1392, 1405,
1409, 1420, 1428, 1439, 1443, 1455, 1468, 1471, 1486, 1487, 1500, 1503, 1516, 1519, 1520, 1523,
1536, 1551, 1557, 1567, 1568, 1583, 1584, 1597, 1603, 1615, 1617, 1629, 1646, 1647, 1648, 1654,
1676, 1679, 1692, 1695, 1708, 1711, 1715, 1725, 1732, 1743, 1744, 1757, 1761, 1770, 1778, 1786,
1804, 1807, 1820, 1823, 1836, 1839, 1854, 1855, 1856, 1868, 1877, 1887, 1902, 1903, 1916, 1919,
1934, 1935, 1942, 1951, 1964, 1967, 1973, 1983, 1988, 1999, 2000, 2013, 2028, 2031, 2034, 2047,
2048, 2051, 2064, 2067, 2080, 2085, 2099, 2104, 2112, 2122, 2131, 2143, 2144, 2147, 2160, 2163,
2177, 2188, 2202, 2207, 2208, 2221, 2227, 2239, 2243, 2254, 2259, 2271, 2279, 2283, 2300, 2303,
2307, 2316, 2320, 2334, 2339, 2348, 2355, 2367, 2369, 2383, 2384, 2395, 2405, 2409, 2416, 2422,
2432, 2435, 2448, 2451, 2464, 2479, 2483, 2493, 2499, 2508, 2512, 2513, 2531, 2537, 2544, 2547,
2563, 2570, 2576, 2591, 2592, 2607, 2611, 2621, 2633, 2636, 2643, 2650, 2656, 2657, 2675, 2679};

inline bool IsUnusedPixel(int elID) {
    return kUnusedCDetPixels.count(elID) != 0;
}

TChain *T = 0;

static const double CDet_y_half_length = 0.30;
static const double ECal_dist = 6.6;

static const double TDC_calib_to_ns = 0.01;

double TDCBinLow;
double TDCBinHigh;
int NTDCBins;
Double_t RefNTDCBins;
Double_t RefLeMin;
Double_t RefLeMax;

std::vector<double> vRefRawLe;
std::vector<double> vRefRawTe;
std::vector<double> vRefRawTot;
std::vector<int>    vRefRawPMT;
std::vector<double> vAllRawLe;
std::vector<double> vAllRawLeNoRef;
std::vector<int>    vAllRawPMT;
std::vector<double> vT_mod;
std::vector<double> vAllRefForLe;

std::vector<std::vector<double>> vCDetPaddleRawTot;
std::vector<std::vector<double>> vCDetPaddleCutTot;

struct CDetHit {
  int    id;     // pixelID
  double le_ns;  // LE
  double tot_ns; // TOT
  double te_ns;  // TE
  double x_pos; //x position in cm
  double y_pos; //y
  double z_pos; //z
};

struct AdjacentHits {
  int    id;     // pixelID
  double le_ns;  // LE
  double tot_ns; // TOT
  double te_ns;  // TE
  double x_pos; //x position in cm
  double y_pos; //y
  double z_pos; //z
};

std::vector<std::vector<CDetHit>> vEventHits; // [event][hit]
std::vector<std::vector<CDetHit>> vAdjEventHits; // [event][hit]

std::vector<int> rawRate(2688, 0); 
int rateEvTrack = 0;
std::vector<double> chanRates(2688,0);
std::vector<int> cutRate(2688, 0); 
int cutRateEvTrack = 0;
std::vector<double> cutChanRates(2688,0);
std::vector<double> ave_tot(2688,0);
std::vector<int> vNumAdjacentHits;
std::vector<double> adjChanRates(2688,0);
std::vector<int> goodRate(2688, 0); 

int NTotBins = 200;
double TotBinLow = 1.;
double TotBinHigh = 51.;

void t_mod(int runnum = 5811, Int_t neventsr=500000, Int_t minSeg = -1, Int_t maxSeg = -1, 
          Double_t LeMin = 0.02, Double_t LeMax = 60, Double_t TotMin = 1, Double_t TotMax = 100,
          Double_t XDiffCut = 0.01, Double_t XOffset = 0.02, Double_t YOffset = 0.1){

    RefLeMin = 0.0;
    RefLeMax = 252.0;
    RefNTDCBins = (RefLeMax-RefLeMin)/4;
    Double_t RefTotMin = 1.0;
    Double_t RefTotMax = 252.0;

    NTDCBins = 2*(LeMax-LeMin)/.0160167; // 4 ns is the trigger time, 0.018 ns is the expected time resolution, if we use a reference TDC ? 
                    // 4 ns resolution is the best we can hope for, I think, using only the module trigger time.
    TDCBinLow = LeMin;
    TDCBinHigh = LeMax;

  if (!T) {
      T = new TChain("T");
      //int onlySegment = -1; // set to >=0 to pick just one
      AddRunFilesToChain(T, REPLAYED_DIR.Data(), runnum, minSeg, maxSeg);
  }

  TTreeReader reader(T);

  //Set TTreeReaders
  /* ----- Earm ----- */ 

  // ******CDet******
  // ----- CDet arrays -----
  TTreeReaderArray<double> TDCmult(reader, "earm.cdet.tdc_mult");

  TTreeReaderArray<double> RawElID   (reader, "earm.cdet.hits.TDCelemID");
  TTreeReaderArray<double> RawElLE   (reader, "earm.cdet.hits.t");
  TTreeReaderArray<double> RawElTE   (reader, "earm.cdet.hits.t_te");
  TTreeReaderArray<double> RawElTot  (reader, "earm.cdet.hits.t_tot");

  TTreeReaderArray<double> GoodElID  (reader, "earm.cdet.hit.pmtnum");
  // TTreeReaderArray<double> GoodElLE  (reader, "earm.cdet.hit.tdc_le");
  // TTreeReaderArray<double> GoodElTE  (reader, "earm.cdet.hit.tdc_te");
  // TTreeReaderArray<double> GoodElTot (reader, "earm.cdet.hit.tdc_tot");

  TTreeReaderArray<double> GoodX     (reader, "earm.cdet.hit.xhit");
  TTreeReaderArray<double> GoodY     (reader, "earm.cdet.hit.yhit");
  TTreeReaderArray<double> GoodZ     (reader, "earm.cdet.hit.zhit");

  // TTreeReaderArray<double> GoodCol   (reader, "earm.cdet.hit.row");
  // TTreeReaderArray<double> GoodRow   (reader, "earm.cdet.hit.col");
  // TTreeReaderArray<double> GoodLayer (reader, "earm.cdet.hit.layer");

  // ----- cdet scalars ----- 
  TTreeReaderValue<double> nhits        (reader, "earm.cdet.nhits");
  // TTreeReaderValue<double> ngoodhits    (reader, "earm.cdet.ngoodhits");
  // TTreeReaderValue<double> ngoodTDChits (reader, "earm.cdet.ngoodTDChits");

  //------ECal-------
  // Cluster arrays
  // TTreeReaderArray<double> ECal_clus_adctime      (reader, "earm.ecal.clus.adctime");
  // TTreeReaderArray<double> ECal_clus_again        (reader, "earm.ecal.clus.again");
  // TTreeReaderArray<double> ECal_clus_atimeblk     (reader, "earm.ecal.clus.atimeblk");
  // TTreeReaderArray<double> ECal_clus_col          (reader, "earm.ecal.clus.col");
  // TTreeReaderArray<double> ECal_clus_e            (reader, "earm.ecal.clus.e");
  // TTreeReaderArray<double> ECal_clus_eblk         (reader, "earm.ecal.clus.eblk");
  // TTreeReaderArray<double> ECal_clus_id           (reader, "earm.ecal.clus.id");
  // TTreeReaderArray<double> ECal_clus_nblk         (reader, "earm.ecal.clus.nblk");
  // TTreeReaderArray<double> ECal_clus_row          (reader, "earm.ecal.clus.row");
  // TTreeReaderArray<double> ECal_clus_x            (reader, "earm.ecal.clus.x");
  // TTreeReaderArray<double> ECal_clus_y            (reader, "earm.ecal.clus.y");

  // Cluster count (scalar)
  //TTreeReaderValue<double> ECal_nclus(reader, "earm.ECal.nclus");

  //event-level ECal branches
  TTreeReaderValue<double> ECalX       (reader, "earm.ecal.x");
  TTreeReaderValue<double> ECalY       (reader, "earm.ecal.y");
  // TTreeReaderValue<double> ECalE       (reader, "earm.ecal.e");
  // TTreeReaderValue<double> ECalAdcTime (reader, "earm.ecal.adctime");

  vCDetPaddleRawTot.assign(2688, std::vector<double>{});
  vCDetPaddleCutTot.assign(2688, std::vector<double>{});

  Int_t Nev = T->GetEntries();
  cout << "N entries in tree is " << Nev << endl;
  Int_t NEventsAnalysis;// = Nev;
  if(neventsr==-1) NEventsAnalysis = Nev;
  else NEventsAnalysis = neventsr;
  cout << "Running analysis for " << NEventsAnalysis << " events" << endl;

  //================================================================= Event Loop
  // variables outside event loop
  Int_t EventCounter = 0;
  cout << "Starting Event Loop" << endl;

  // event loop start
  Int_t event = 0;
  while(reader.Next()){
    event++;
    event = event - 1;
    EventCounter++;
    // Only stop early if neventsr > 0
    if (neventsr > 0 && EventCounter > neventsr) {
        break;
    }

    Int_t nh = *nhits;
    
    if (EventCounter % 1000 == 0) {
      cout << EventCounter << "/" << NEventsAnalysis << "/ Nhits = " << (Int_t)nh << endl;
    }

    // First pass through hits:  purpose is to get reference LE TDC Value for this event
    
    double event_ref_tdc = 0.0;
    double ref_int = 0;
    double ref_corr = 0;
    // std::cout << "RawElID Size = " << " " << RawElID.GetSize() << std::endl;
    // for (auto val : RawElLE) {
    //   std::cout << " rawLE array= "<< " " << val <<std::endl; 
    // }
    for(Int_t el=0; el<RawElID.GetSize(); el++) {
      if ((Int_t)RawElID[el] == 2696) {  // only look at ref PMT 
        bool good_ref_le_time = RawElLE[el] > 0.0/TDC_calib_to_ns && RawElLE[el] <= 252.0/TDC_calib_to_ns;
        bool good_ref_event = good_ref_le_time;
        if ( good_ref_event ) {
          //std::cout << "event = " << " " << EventCounter << " " << "ref time = " << RawElLE[el] << std::endl;
          //std::cout << "ref Tot = " << RawElTot[el] << std::endl;
          if ( (Int_t)RawElID[el] == 2696 && (Int_t)RawElLE[el] > 0 ) {
              vRefRawLe.push_back(RawElLE[el] * TDC_calib_to_ns);
              vRefRawTe.push_back(RawElTE[el] * TDC_calib_to_ns);
              vRefRawTot.push_back(RawElTot[el] * TDC_calib_to_ns);
              vRefRawPMT.push_back((int)RawElID[el]);

              event_ref_tdc = RawElLE[el]*TDC_calib_to_ns - 48;
              ref_int = std::floor(event_ref_tdc);
              ref_corr = event_ref_tdc - ref_int;
          }
        }
      }
    }// end ref TDC loop
      
    rateEvTrack++;
    std::vector<CDetHit> eventHits;
    int CDetPassedBoolCount = 0;
    // Build lookup from PMT id -> index in Good* arrays for this TTree entry
    std::unordered_map<int,int> goodIdx;
    goodIdx.reserve(GoodElID.GetSize());
    for (int ig = 0; ig < (int)GoodElID.GetSize(); ig++) {
      goodIdx[(int)GoodElID[ig]] = ig;
    }

    for(Int_t el=0; el<RawElID.GetSize(); el++){

      const int raw_pmt = (int)RawElID[el];
      auto itGood = goodIdx.find(raw_pmt);
      const bool hasGood = (itGood != goodIdx.end());
      const int ig = hasGood ? itGood->second : -1;
      const double gx = hasGood ? GoodX[ig] : 1.0e9;
      const double gy = hasGood ? GoodY[ig] : 1.0e9;
      const double gz = hasGood ? GoodZ[ig] : 1.0e9;

      if (kUnusedCDetPixels.count(RawElID[el])) continue; //checks if pixel is unused, if it is, skip. 
      bool good_raw_le_time = RawElLE[el] >= LeMin/TDC_calib_to_ns && RawElLE[el] <= LeMax/TDC_calib_to_ns;
      bool goodhit_tot = RawElTot[el] >= TotMin/TDC_calib_to_ns && RawElTot[el] <= TotMax/TDC_calib_to_ns;
      // bool good_ECal_diff_x = (gx-((*ECalX)*(gz)/ECal_dist)-XOffset) <= XDiffCut && 
      //     (gx-((*ECalX)*(gz)/ECal_dist)-XOffset) >= -1.0*XDiffCut;
      // bool good_ECal_diff_y = (gy-((*ECalY)*(gz)/ECal_dist)-YOffset) <= 1.2*CDet_y_half_length && 
      //     (gy-((*ECalY)*(gz)/ECal_dist)-YOffset) >= -1.2*CDet_y_half_length;
      bool good_raw_event = good_raw_le_time && goodhit_tot;// && good_ECal_diff_x && good_ECal_diff_y;
      

    //if ((Int_t)RawElID[el] > 1000) cout << "el = " << el << " Hit ID = " << (Int_t)RawElID[el] << "    TDC = " << RawElLE[el]*TDC_calib_to_ns << endl;
    //cout << "Raw ID = " << RawElID[el] << " raw le = " << RawElLE[el] << " raw te = " << RawElTE[el] << " raw tot = " << RawElTot[el] << endl;
      if ( good_raw_event ) {
        CDetPassedBoolCount++;
        int idx = RawElID[el];
        if (0 <= idx && idx < 2688) {
          double tot_ns = RawElTot[el]*TDC_calib_to_ns;
          double le_ns = RawElLE[el]*TDC_calib_to_ns - event_ref_tdc;
          double te_ns = RawElTE[el]*TDC_calib_to_ns - event_ref_tdc;
          rawRate[idx]++;
          vCDetPaddleRawTot[idx].push_back(tot_ns);
          eventHits.push_back({idx, le_ns, tot_ns, te_ns, gx, gy, gz});
        } //getting rates and tot for pixels
        if ( (Int_t)RawElID[el] < 2688 ) {
          //fill all hits vectors
          double t_ns = RawElLE[el]*TDC_calib_to_ns - event_ref_tdc;//ref_corr;
          double wrapped = fmod(t_ns,2.0);

          vAllRawLe.push_back(RawElLE[el]*TDC_calib_to_ns - event_ref_tdc);//ref_corr);
          vAllRawLeNoRef.push_back(RawElLE[el]*TDC_calib_to_ns);//ref_corr);
          vAllRawPMT.push_back(RawElID[el]);
          vT_mod.push_back(wrapped);
          vAllRefForLe.push_back(event_ref_tdc);
        }
      }
    }// all raw tdc hit loop
    if (CDetPassedBoolCount >= 1){
      vEventHits.push_back(eventHits);
    }
    //check nadjacent pairs for each event (pixels 260-270)
    auto is_unused = [&](int id) -> bool {
      return kUnusedCDetPixels.count(id) != 0;
    };

    auto is_adjacent_with_skip = [&](int a, int b) -> bool {
      if (b == a + 1) {
        // adjacent normally, but only if that neighbor isn't unused
        return !is_unused(b);
      }
      if (b == a + 2) {
        // treat as adjacent if the in-between pixel is unused
        return is_unused(a + 1) && !is_unused(b);
      }
      return false;
    };

    if (CDetPassedBoolCount >= 1) {
      std::unordered_set<int> adjIDs; //per-event
      std::vector<int> ids;
      const auto& currentEvent = vEventHits.back();
      ids.reserve(currentEvent.size());
      adjIDs.reserve(currentEvent.size());

      for (const auto& hit : currentEvent) {
        int id = hit.id;
        // if (260 <= id && id <= 270) ids.push_back(id);
        ids.push_back(id);
      }

      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

      int nAdjacentHits = 0;
      for (size_t i = 0; i + 1 < ids.size(); i++) {
        int a = ids[i];
        int b = ids[i+1];
        if (is_adjacent_with_skip(a,b)) {
          adjIDs.insert(a);
          adjIDs.insert(b);
          nAdjacentHits++;
        }
      }

      for (int id : adjIDs) {
        if (0 <= id && id < (int)goodRate.size()){
          goodRate[id]++;
        }
      }

      vNumAdjacentHits.push_back(nAdjacentHits);

      //filer hit list for adj pairs only
      std::vector<CDetHit> adjHitsThisEvent;
      adjHitsThisEvent.reserve(currentEvent.size());

      for (const auto& hit : currentEvent){
        if (adjIDs.count(hit.id)){
          adjHitsThisEvent.push_back(hit);
        }
      }
      vAdjEventHits.push_back(std::move(adjHitsThisEvent)); //should contain adj pairs only
    }

    //std::cout << "Event = " << rateEvTrack << std::endl;
    if (CDetPassedBoolCount >= 1){
      const auto& currEvent = vEventHits.back();
      for (const auto& hit : currEvent){
        if (rateEvTrack <= 5){
          int id = hit.id;
          double le = hit.le_ns;
          double tot = hit.tot_ns;
          //std::cout << " Pixel ID = " << id << ", LE = " << le << ", TOT = " << tot << std::endl;
        }
      }
    }

  }//end event loop 
  // std::cout << "nevents = " << rateEvTrack << std::endl;
  for (int i = 0; i < 2688; i++){
    chanRates[i] = (double)rawRate[i] / (rateEvTrack);
    adjChanRates[i] = (double)goodRate[i] / (rateEvTrack);
    //std::cout << "triggered Rate in Pixel " << 417 + i << " = " << chanRates << " & with time window Rate = " << chanRates / winWidth <<std::endl;
  }

  //Second Pass over all events for tot_ave calc

  for (Int_t idx = 0; idx < 2688; idx++){
    Int_t ihits = vCDetPaddleRawTot[idx].size();
    double sumTot = 0;
    for (Int_t i = 0; i < ihits; i++){
      sumTot += vCDetPaddleRawTot[idx][i];
    }
    ave_tot[idx] = sumTot / ihits;
  }

  for (Int_t idx = 0; idx < 2688; idx++){
    Int_t ihits = vCDetPaddleRawTot[idx].size();
    for (Int_t i = 0; i < ihits; i++){
      double hit_tot = vCDetPaddleRawTot[idx][i];
      if (hit_tot > ave_tot[idx]){
        vCDetPaddleCutTot[idx].push_back(hit_tot);
        cutRate[idx]++;
      }
    }
    cutChanRates[idx] = (double)cutRate[idx] / rateEvTrack;
  }
}//end main

void plotAdjTiming(double width = 4, double leMin=0, double leMax = 65, double totMin = 0, double totMax = 150){
  int TDCBinNum = (int)((leMax-leMin)/width);
  TH1D* hLe = new TH1D("hLe", "Leading Edge Pixel; LE (ns);Counts",TDCBinNum,leMin,leMax);
  TH1D* hTe = new TH1D("hTe", "Trailing Edge Pixel; TE (ns);Counts",TDCBinNum,leMin,leMax+totMax);
  TH1D* hTot = new TH1D("hTot", "Tot Pixel; Tot (ns);Counts",TDCBinNum,totMin,totMax);
  TH2D* hLEvsTOT = new TH2D("hLEvsTOT", "LE vs TOT Pixel; TOT (ns);LE (ns)",TDCBinNum,totMin,totMax,TDCBinNum,leMin,leMax); 
  TH2D* hTEvsLE = new TH2D("hTEvsLE", "TE vs LE Pixel; LE (ns);TE (ns)",TDCBinNum,leMin,leMax,TDCBinNum,0,70);

  for (const auto& event : vAdjEventHits){
    for (const auto& hit : event){
      double tot = hit.tot_ns;
      double le = hit.le_ns;
      double te = hit.te_ns;
      hLe->Fill(le);
      hTe->Fill(te);
      hTot->Fill(tot);
      hLEvsTOT->Fill(tot,le);
      hTEvsLE->Fill(le,te);
    }
  }

  //define canvas and plot
  TCanvas* cTiming = new TCanvas("cTiming", "Timing Plots", 900,700);
  cTiming->Divide(1,3); //2x3
  //LE
  cTiming->cd(1);
  hLe->Draw();
  //TE
  cTiming->cd(2);
  hTe->Draw();
  //TOT
  cTiming->cd(3);
  hTot->Draw();

  TCanvas* cLEvsTOT = new TCanvas("cLEvsTOT", "Leading Edge vs Tot",900,700);
  hLEvsTOT->Draw("COLZ");

  TCanvas* cTEvsLE = new TCanvas("cTEvsLE", "TE vs LE",900,700);
  hTEvsLE->Draw("COLZ");

}

void plotLowRate(double rate = 0.3, bool yesTotCut = true, double Width = 4, double leMin = 0, double leMax = 60, double totMin = 0, double totMax = 100){
  TH1::AddDirectory(kFALSE);
  int TDCBinNum = (int)((leMax-leMin)/Width);
  TH1D* hLe1 = new TH1D("hLe1", "Leading Edge Pixel; LE (ns);Counts",TDCBinNum,leMin,leMax);
  TH1D* hTe1 = new TH1D("hTe1", "Trailing Edge Pixel; TE (ns);Counts",TDCBinNum,leMin,leMax+totMax);
  TH1D* hTot1 = new TH1D("hTot1", "Tot Pixel; Tot (ns);Counts",TDCBinNum,totMin,totMax);
  TH2D* hLEvsTOT1 = new TH2D("hLEvsTOT1", "LE vs TOT Pixel; TOT (ns);LE (ns)",TDCBinNum,totMin,totMax,TDCBinNum,leMin,leMax); 
  TH2D* hTEvsLE1 = new TH2D("hTEvsLE1", "TE vs LE Pixel; LE (ns);TE (ns)",TDCBinNum,leMin,leMax,TDCBinNum,0,70);
  
  //loop through hits and events
  for (const auto& event : vEventHits){
    for (const CDetHit& hit : event){
      int id = hit.id;
      double totCut = 0;
      if (yesTotCut) totCut = ave_tot[id];
      if (chanRates[id] <= rate && hit.tot_ns > totCut){  
        double tot = hit.tot_ns;
        double le = hit.le_ns;
        double te = hit.te_ns;
        hLe1->Fill(le);
        hTe1->Fill(te);
        hTot1->Fill(tot);
        hLEvsTOT1->Fill(tot,le);
        hTEvsLE1->Fill(le,te);
      }
    }//loop through hits
  }//loop through events

  //define canvas and plot
  TCanvas* cTiming = new TCanvas("cTiming", "Timing Plots", 900,700);
  cTiming->Divide(1,3); //2x3
  //LE
  cTiming->cd(1);
  hLe1->Draw();
  //TE
  cTiming->cd(2);
  hTe1->Draw();
  //TOT
  cTiming->cd(3);
  hTot1->Draw();

  TCanvas* cLEvsTOT = new TCanvas("cLEvsTOT", "Leading Edge vs Tot",900,700);
  hLEvsTOT1->Draw("COLZ");

  TCanvas* cTEvsLE = new TCanvas("cTEvsLE", "TE vs LE",900,700);
  hTEvsLE1->Draw("COLZ");

}

void plotNumAdjacent(int nbins = 50){
  TH1::AddDirectory(kFALSE);
  TH1D* hNumAdjacentHits = new TH1D("hNumAdjacentHits", "Number Hits in Adjacent Pixels", nbins, 0, nbins);

  for (const auto& hit : vNumAdjacentHits){
    if (hit > 0) hNumAdjacentHits->Fill(hit);
  }

  TCanvas* cNumAdjacentHits = new TCanvas("cNumAdjacentHits", "Number of Adjacent Hits", 900,700);
  hNumAdjacentHits->Draw();

}

void plotAveTotPerPixel() {
  const int nPixels = ave_tot.size();

  TH1::AddDirectory(kFALSE);

  TH1D* hAveTot = new TH1D("hAveTot","Average TOT per CDet pixel;Pixel ID;Average TOT (ns)",nPixels, -0.5, nPixels - 0.5);

  for (int p = 0; p < nPixels; ++p) {
    if (ave_tot[p] > 0) {   // optional guard
      hAveTot->SetBinContent(p + 1, ave_tot[p]);
    }
  }

  TCanvas* cAveTot = new TCanvas("cAveTot", "Average TOT per Pixel", 1200, 500);
  hAveTot->Draw("HIST");
}

void plotPixelComp(int pixel1 = 257, int pixel2 = 261, double Width = 1, double leMin = 0, double leMax = 60, double totMin = 0, double totMax = 100){
  TH1::AddDirectory(kFALSE);
  int TDCBinNum = (int)((leMax-leMin)/Width);
  TH1D* hLe1 = new TH1D("hLe1", Form("Leading Edge Pixel %d; LE (ns);Counts", pixel1),TDCBinNum,leMin,leMax);
  TH1D* hTe1 = new TH1D("hTe1", Form("Trailing Edge Pixel %d; TE (ns);Counts", pixel1),TDCBinNum,leMin,leMax+totMax);
  TH1D* hTot1 = new TH1D("hTot1", Form("Tot Pixel %d; Tot (ns);Counts", pixel1),TDCBinNum,totMin,totMax);
  TH1D* hLe2 = new TH1D("hLe2", Form("Leading Edge Pixel %d; LE (ns);Counts", pixel2),TDCBinNum,leMin,leMax);
  TH1D* hTe2 = new TH1D("hTe2", Form("Trailing Edge Pixel %d; TE (ns);Counts", pixel2),TDCBinNum,leMin,leMax+totMax);
  TH1D* hTot2 = new TH1D("hTot2", Form("Tot Pixel %d; Tot (ns);Counts", pixel2),TDCBinNum,totMin,totMax);
  TH2D* hLEvsTOT1 = new TH2D("hLEvsTOT1", Form("LE vs TOT Pixel %d; TOT (ns);LE (ns)", pixel1),TDCBinNum,totMin,totMax,TDCBinNum,leMin,leMax);
  TH2D* hLEvsTOT2 = new TH2D("hLEvsTOT2", Form("LE vs TOT Pixel %d; TOT (ns);LE (ns)", pixel2),TDCBinNum,totMin,totMax,TDCBinNum,leMin,leMax);
  TH2D* hTEvsLE1 = new TH2D("hTEvsLE1", Form("TE vs LE Pixel %d; TE (ns);LE (ns)", pixel1),TDCBinNum,leMin,leMax,TDCBinNum,0,70);
  TH2D* hTEvsLE2 = new TH2D("hTEvsLE2", Form("TE vs LE Pixel %d; TE (ns);LE (ns)", pixel2),TDCBinNum,leMin,leMax,TDCBinNum,0,70);

  //pick out pixels
  for (const auto& event : vEventHits){
    for (const CDetHit& hit : event){
      int id = hit.id;
      if (id == pixel1 && hit.tot_ns > totMin){
        double tot = hit.tot_ns;
        double le = hit.le_ns;
        double te = hit.te_ns;
        hLe1->Fill(le);
        hTe1->Fill(te);
        hTot1->Fill(tot);
        hLEvsTOT1->Fill(tot,le);
        hTEvsLE1->Fill(le,te);
      }
      else if (id == pixel2 && hit.tot_ns > totMin){
        double tot = hit.tot_ns;
        double le = hit.le_ns;
        double te = hit.te_ns;
        hLe2->Fill(le);
        hTe2->Fill(te);
        hTot2->Fill(tot);
        hLEvsTOT2->Fill(tot,le);
        hTEvsLE2->Fill(le,te);
      }
    }//loop throuh hits
  }//loop through events

  //define canvas and plot
  TCanvas* cTiming = new TCanvas("cTiming", "Timing Plots", 900,700);
  cTiming->Divide(2,3); //2x3
  //LE
  cTiming->cd(1);
  hLe1->Draw();
  cTiming->cd(2);
  hLe2->Draw();
  //TE
  cTiming->cd(3);
  hTe1->Draw();
  cTiming->cd(4);
  hTe2->Draw();
  //TOT
  cTiming->cd(5);
  hTot1->Draw();
  cTiming->cd(6);
  hTot2->Draw();

  TCanvas* cLEvsTOT = new TCanvas("cLEvsTOT", "Leading Edge vs Tot",900,700);
  cLEvsTOT->Divide(1,2);
  cLEvsTOT->cd(1);
  hLEvsTOT1->Draw();
  cLEvsTOT->cd(2);
  hLEvsTOT2->Draw();

  TCanvas* cTEvsLE = new TCanvas("cTEvsLE", "TE vs LE",900,700);
  cTEvsLE->Divide(1,2);
  cTEvsLE->cd(1);
  hTEvsLE1->Draw();
  cTEvsLE->cd(2);
  hTEvsLE2->Draw();
}

void plotSingleTot(int pixel_base = 0, bool raw = true, double width = 1, double totMin=1, double totMax=80){
  TH1::AddDirectory(kFALSE);
  if (pixel_base % 16 != 0) {
    Error("plotSingleTot", "pixel_base = %d is not a multiple of 16", pixel_base);
    return;
  }

  const int nPlots = 16;
  int TDCBinNum = (int)((totMax-totMin)/width);

  // Decode pixel → {layer, side, submodule, pmt, pixel}
  auto info = getLocation(pixel_base);

  int layer     = info[0] + 1;  // display as 1-based
  int side      = info[1];      // 0=L, 1=R
  int submodule = info[2] + 1;
  int bar       = info[3] + 1;

  TString sideStr = (side == 0) ? "L" : "R";

  TString canvasTitle = Form("Layer %d | %s | Module %d | Bar %d", layer, sideStr.Data(), submodule, bar);
  // Canvas with 4x4 pads
  TString cname = Form("cTot_%d", pixel_base);
  TCanvas* cTot = new TCanvas(cname, canvasTitle, 1200, 1000);
  cTot->Divide(4, 4, 0.001, 0.001);

  // Histogram array
  TH1D* hTot[nPlots];

  for (int i = 0; i < nPlots; i++) {

    int pixel = pixel_base + i;

    TString hname  = Form("hTot_pix%d", pixel);
    if (raw){
      TString htitle = Form("Pixel %d;TOT (ns);Counts", pixel);
      hTot[i] = new TH1D(hname, htitle, TDCBinNum, totMin, totMax);
    }
    if (!raw){
      TString htitle = Form("Pixel %d w/ TOT > %f;TOT (ns);Counts", pixel, ave_tot[pixel]);
      hTot[i] = new TH1D(hname, htitle, TDCBinNum, totMin, totMax);
    }

    // Fill histogram
    if (raw){
      for (const auto& x : vCDetPaddleRawTot[pixel]) {
        hTot[i]->Fill(x);
      }
    }
    if (!raw){
      for (const auto& x : vCDetPaddleCutTot[pixel]) {
        hTot[i]->Fill(x);
      }
    }

    // Draw
    cTot->cd(i + 1);
    hTot[i]->Draw();

    // If unused pixel: draw a black square in the top-right corner of the pad
    if (kUnusedCDetPixels.count(pixel)) {
      // NDC coordinates: (x1,y1,x2,y2) in [0,1] pad coordinates
      TPaveText* flag = new TPaveText(0.82, 0.82, 0.95, 0.95, "NDC");
      flag->SetFillColor(kBlack);
      flag->SetLineColor(kBlack);
      flag->SetBorderSize(1);
      flag->AddText("");       // empty; just a filled box
      flag->Draw("same");
    }
  }

  cTot->Update();
}

void getRate(int pixel, bool cut = false){
  if (!cut) std:: cout << "Rate in Pixel " << pixel << " = " << chanRates[pixel] << std::endl;
  if (cut) std:: cout << "Rate in Pixel " << pixel << " (with Tot Cut) = " << cutChanRates[pixel] << std::endl;
}

void plotRateVsID(bool raw = true){
  TH1::AddDirectory(kFALSE);
    // --- constants ---
  const int NCHAN_TOTAL = 2688;
  const int NCHAN_LAYER = 1344;
  const int NCHAN_SIDE  = 672;
  const int NMOD        = 3;
  const int NCHAN_MOD   = NCHAN_SIDE / NMOD; // 224

  // Helper: create one segment histogram and fill from chanRates
  auto MakeRateHist = [&](const char* hname,
                          const char* htitle,
                          int idStart, int idEnd,
                          double yMax = 1.5) -> TH1D* {
    const int nbins = idEnd - idStart + 1; // inclusive
    TH1D* h = new TH1D(hname, htitle, nbins, idStart, idEnd + 1); // [start, end+1)
    h->SetStats(0);
    h->SetMinimum(0.0);
    h->SetMaximum(yMax);

    for (int id = idStart; id <= idEnd; id++) {
      const int bin = h->FindBin(id);
      if (raw){
        h->SetBinContent(bin, chanRates[id]);
      }
      if (!raw){
        h->SetBinContent(bin, adjChanRates[id]);
      }
    }
    return h;
  };

  struct Seg { int layer; int mod; const char* side; int start; int end; };

  std::vector<Seg> segs;
  segs.reserve(12);

  auto AddLayerSegs_LeftThenRight = [&](int layer, int base) {
    const int L0 = base + 0;
    const int R0 = base + NCHAN_SIDE;

    // Left side: M1, M2, M3
    for (int m = 0; m < NMOD; m++) {
      int s = L0 + m*NCHAN_MOD;
      int e = s + NCHAN_MOD - 1;
      segs.push_back({layer, m+1, "L", s, e});
    }
    // Right side: M1, M2, M3
    for (int m = 0; m < NMOD; m++) {
      int s = R0 + m*NCHAN_MOD;
      int e = s + NCHAN_MOD - 1;
      segs.push_back({layer, m+1, "R", s, e});
    }
  };

  // Layer 1: IDs 0..1343
  AddLayerSegs_LeftThenRight(1, 0);
  // Layer 2: IDs 1344..2687
  AddLayerSegs_LeftThenRight(2, 1344);

  // --- build histograms ---
  TH1D* hRateSeg[12] = {nullptr};

  for (int i = 0; i < 12; i++) {
    const auto& s = segs[i];
    if (raw) {
      TString name  = Form("hRateVsIDL%dM%d%s", s.layer, s.mod, s.side);
      TString title = Form("CDet L%d %s M%d Rate vs Pixel ID;Pixel ID;Rate",
                          s.layer, s.side, s.mod);
      hRateSeg[i] = MakeRateHist(name.Data(), title.Data(), s.start, s.end, 1.5);
    }
    if (!raw){
      TString name  = Form("hRateVsIDL%dM%d%s", s.layer, s.mod, s.side);
      TString title = Form("CDet L%d %s M%d Rate w/Adj vs Pixel ID;Pixel ID;Rate",
                          s.layer, s.side, s.mod);
      hRateSeg[i] = MakeRateHist(name.Data(), title.Data(), s.start, s.end, 1.5);
    }
  }

  // --- Draw: Layer 1 canvas (Left M1-3 then Right M1-3) ---
  TCanvas* cRateL1 = new TCanvas("cRateL1", "CDet Rate vs ID (Layer 1)", 1400, 800);
  cRateL1->Divide(3,2); // top row: left M1-3, bottom row: right M1-3

  int pad = 1;
  for (int i = 0; i < 12; i++) {
    if (segs[i].layer != 1) continue;
    cRateL1->cd(pad++);
    hRateSeg[i]->Draw("HIST");

    int nb = hRateSeg[i]->GetNbinsX();
    double avg = hRateSeg[i]->Integral(1,nb) / nb;
    double xmin = hRateSeg[i]->GetXaxis()->GetXmin();
    double xmax = hRateSeg[i]->GetXaxis()->GetXmax();

    TLine* lAvg = new TLine(xmin, avg, xmax, avg);
    lAvg->SetLineColor(kRed);
    lAvg->SetLineStyle(2);
    lAvg->SetLineWidth(2);
    lAvg->Draw("SAME");
  }

  // --- Draw: Layer 2 canvas (Left M1-3 then Right M1-3) ---
  TCanvas* cRateL2 = new TCanvas("cRateL2", "CDet Rate vs ID (Layer 2)", 1400, 800);
  cRateL2->Divide(3,2);

  pad = 1;
  for (int i = 0; i < 12; i++) {
    if (segs[i].layer != 2) continue;
    cRateL2->cd(pad++);
    hRateSeg[i]->Draw("HIST");

    int nb = hRateSeg[i]->GetNbinsX();
    double avg = hRateSeg[i]->Integral(1,nb) / nb;
    double xmin = hRateSeg[i]->GetXaxis()->GetXmin();
    double xmax = hRateSeg[i]->GetXaxis()->GetXmax();

    TLine* lAvg = new TLine(xmin, avg, xmax, avg);
    lAvg->SetLineColor(kRed);
    lAvg->SetLineStyle(2);
    lAvg->SetLineWidth(2);
    lAvg->Draw("SAME");
  }
}

TCanvas *plotAllTDC(double width = 4, double TDCBinLow = 0, double TDCBinHigh = 65){
    //define histograms
    int TDCBinNum = (int)((TDCBinHigh-TDCBinLow)/width);
    TH1F *hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
    TH1F *hT_mod = new TH1F("hT_mod","hT_mod",
            200, 0, 2.5);
    //fill necessary histograms from vectors
    for (double x : vAllRawLe) hAllRawLe->Fill(x);
    for (double x : vT_mod) hT_mod->Fill(x);

    TCanvas *c = new TCanvas("c", "Tmod (0,2)",800,800);
    // c->Divide(2,1);

    c->cd(1);
    hAllRawLe->Draw();

    // c->cd(2);
    // hT_mod->Draw("Fill");
    return c;
}

TCanvas *plotNoRef(double TDCBinLow, double TDCBinHigh){
    //define histograms
    TH1F *hAllRawLeNoRef = new TH1F("hAllRawLeNoRef",
            "hAllRawLeNoRef",
            NTDCBins, TDCBinLow, TDCBinHigh);
    //fill necessary histograms from vectors
    for (double x : vAllRawLeNoRef) hAllRawLeNoRef->Fill(x);

    TCanvas *c = new TCanvas("c", "LE no Ref",800,800);
    hAllRawLeNoRef->Draw();
    return c;
}

TCanvas *plotRef(Double_t cutTotMin = 90, Double_t cutTotMax = 150, Int_t RefNTDCBins = 252, Double_t RefLeMin = 0.0, Double_t RefLeMax = 52.0, Int_t RefNTotBins= 252, Double_t RefTotMin = 1, Double_t RefTotMax = 250 ){
    //define histograms
    TH1F *hRefLe = new TH1F("hRefLe",
            "hRefLe",
            RefNTDCBins, RefLeMin, RefLeMax);
    TH1F *hRefTe = new TH1F("hRefTe",
            "hRefTe",
            RefNTDCBins, RefLeMin+RefTotMin, RefLeMax+RefTotMax);
    TH1F *hRefTot = new TH1F("hRefTot",
            "hRefTot",
            RefNTotBins, RefTotMin, RefTotMax);
  
    size_t N = vRefRawLe.size();
    for (size_t i = 0; i < N; i++){
      if (vRefRawTot[i] >= cutTotMin && vRefRawTot[i] <= cutTotMax) {
        hRefLe->Fill(vRefRawLe[i]);
        hRefTe->Fill(vRefRawTe[i]);
        hRefTot->Fill(vRefRawTot[i]);
      }
    }

    TCanvas *c = new TCanvas("c", "Ref Signal",800,800);
    c->Divide(1,3);

    c->cd(1);
    //gPad->SetLogy();
    hRefLe->Draw("HIST");

    c->cd(2);
    //gPad->SetLogy();
    hRefTe->Draw("HIST");

    c->cd(3);
    //gPad->SetLogy();
    hRefTot->Draw("HIST");

    return c;
}

TCanvas *plotRefvsLE(){
  TH2F *hRefVsLE = new TH2F("hRefVsLE", "hRefvsLE;LE (ns);Ref time (ns)",NTDCBins,TDCBinLow,TDCBinHigh,RefNTDCBins, RefLeMin, RefLeMax);

  for (size_t i = 0; i < vAllRawLe.size(); ++i){
    hRefVsLE->Fill(vAllRawLe[i],vAllRefForLe[i]);
  }
  TCanvas* c1 = new TCanvas("c2", "Ref vs LE", 800,800);
  hRefVsLE->Draw("COLZ");
  return c1;
}
