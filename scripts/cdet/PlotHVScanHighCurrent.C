#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>
#include <fstream>
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

std::vector<TCanvas*> canvas_vector;

static const int TDCmult_cut = 100;
static const double xcut = 998.0;
//static const int nhitcutlow = 2;
//static const int nhitcuthigh = 20;
static const double TDC_calib_to_ns = 0.01;
static const double HotChannelRatio = .01;

static const int NumPaddles = 16;
static const int NumBars = 14;
static const int NumLayers = 2;
static const int NumSides = 2;
static const int NumModules = 3;
static const int NumHalfModules = NumModules*NumSides*NumLayers;

static const int NumCDetPaddles = NumHalfModules*NumBars*NumPaddles; //2688
static const int nRef = 4;
static const int NumRefPaddles = 4;
static const int nTdc = NumCDetPaddles+NumRefPaddles; //2704

static const int NumSidesTotal = NumSides*NumLayers;
static const int NumCDetPaddlesPerSide = NumCDetPaddles/NumSidesTotal; //672
static const int NumLogicalPaddlesPerSide = NumCDetPaddlesPerSide+nRef; //676

static const int nBarsADC = 0;
static const double ADCCUT = 150.;   //100.0

static const double ecal_dist = 6.6;
static const double cdet_dist_offset = 2.0;
static const double cdet_y_half_length = 0.30;

// List of x-positions (or bins) for unused pixels
static std::vector<double> missingPixelBins = {3, 13, 28, 31, 41, 42, 57, 59, 65, 79, 83, 95, 109, 111, 115, 127, 140, 143, 145, 
  156, 172, 175, 176, 188, 195, 199, 213, 220, 236, 239, 244, 255, 268, 271, 284, 287, 300, 303, 307, 319, 332, 335, 339, 
  351, 354, 364, 371, 381, 384, 396, 401, 410, 419, 423, 435, 436, 451, 461, 465, 479, 480, 483, 508, 511, 512, 515, 540, 
  543, 546, 559, 563, 573, 576, 589, 596, 605, 609, 610, 627, 638, 643, 655, 656, 665, 675, 674, 703, 696, 709, 707, 729, 
  725, 748, 738, 766, 752, 780, 777, 791, 784, 812, 800, 828, 818, 847, 844, 860, 850, 868, 867, 885, 884, 904, 900, 927, 
  912, 943, 940, 947, 945, 971, 967, 991, 986, 1007, 1005, 1023, 1011, 1028, 1027, 1050, 1043, 1068, 1066, 1075, 1072, 1102, 
  1088, 1119, 1106, 1135, 1121, 1151, 1148, 1166, 1162, 1182, 1178, 1186, 1184, 1215, 1203, 1231, 1228, 1247, 1235, 1252, 
  1249, 1274, 1267, 1285, 1282, 1310, 1299, 1321, 1317, 1341, 1340, 1349, 1359, 1372, 1375, 1376, 1391, 1392, 1405, 1409, 
  1420, 1428, 1439, 1443, 1455, 1468, 1471, 1486, 1487, 1500, 1503, 1516, 1519, 1520, 1523, 1536, 1551, 1557, 1567, 1568, 
  1583, 1584, 1597, 1603, 1615, 1617, 1629, 1646, 1647, 1648, 1654, 1676, 1679, 1692, 1695, 1708, 1711, 1715, 1725, 1732, 
  1743, 1744, 1757, 1761, 1770, 1778, 1786, 1804, 1807, 1820, 1823, 1836, 1839, 1854, 1855, 1856, 1868, 1877, 1887, 1902, 
  1903, 1916, 1919, 1934, 1935, 1942, 1951, 1964, 1967, 1973, 1983, 1988, 1999, 2000, 2013, 2031, 2028, 2047, 2034, 2051, 
  2048, 2067, 2064, 2085, 2080, 2104, 2099, 2122, 2112, 2143, 2131, 2147, 2144, 2163, 2160, 2188, 2177, 2207, 2202, 2221, 
  2208, 2239, 2227, 2254, 2243, 2271, 2259, 2283, 2279, 2303, 2300, 2316, 2307, 2334, 2320, 2348, 2339, 2367, 2355, 2383, 
  2369, 2395, 2384, 2409, 2405, 2422, 2416, 2435, 2432, 2451, 2448, 2479, 2464, 2493, 2483, 2508, 2499, 2513, 2512, 2537, 
  2531, 2547, 2544, 2570, 2563, 2591, 2576, 2607, 2592, 2621, 2611, 2636, 2633, 2650, 2643, 2657, 2656, 2679, 2675};

const TString REPLAYED_DIR = "/adaqfs/home/a-onl/sbs/Rootfiles";
const TString ANALYSED_DIR = "/adaqfs/home/a-onl/sbs/Rootfiles/cdetFiles/cdet_histfiles";

// // for local analysis at uog (please leave in comments)
// TString REPLAYED_DIR = "/w/work0/home/rachel/HallA/BB_Hodo/FallRun2021/Replayed";
//TString REPLAYED_DIR = "/w/work2/jlab/halla/sbs_hodo/Rootfiles";
//TString ANALYSED_DIR = "/w/work2/jlab/halla/sbs_hodo/Rootfiles/bbhodo_hist";
//TString ANALYSED_DIR = "/w/work0/home/rachel/HallA/BB_Hodo/FallRun2021/Analysed";

namespace TCDet {
  Int_t NdataMult;
  Double_t TDCmult[nTdc*2];

  Int_t NdataRawElID;
  Double_t RawElID[nTdc*2];
  Int_t NdataRawElLE;
  Double_t RawElLE[nTdc*2];
  Int_t NdataRawElTE;
  Double_t RawElTE[nTdc*2];
  Int_t NdataRawElTot;
  Double_t RawElTot[nTdc*2];
  
  Int_t NdataGoodRow;
  Double_t GoodRow[nTdc*2];
  Int_t NdataGoodCol;
  Double_t GoodCol[nTdc*2];
  Int_t NdataGoodLayer;
  Double_t GoodLayer[nTdc*2];

  Int_t NdataGoodElID;
  Double_t GoodElID[nTdc*2];
  Int_t NdataGoodElLE;
  Double_t GoodElLE[nTdc*2];
  Int_t NdataGoodElTE;
  Double_t GoodElTE[nTdc*2];
  Int_t NdataGoodElTot;
  Double_t GoodElTot[nTdc*2];

  Int_t NdataGoodX;
  Double_t GoodX[nTdc*2];
  Int_t NdataGoodY;
  Double_t GoodY[nTdc*2];
  Int_t NdataGoodZ;
  Double_t GoodZ[nTdc*2];

  Double_t GoodECalX;
  Double_t GoodECalY;
  Double_t GoodECalE;
  Double_t nhits;
  Double_t ngoodhits;
  Double_t ngoodTDChits;

  Double_t nhits_paddles[nTdc*2];
  Double_t ngoodhits_paddles[nTdc*2];
  Double_t ngoodTDChits_paddles[nTdc*2];
  Double_t npaddles;
  Double_t ngoodpaddles;
  Double_t ngoodTDCpaddles;

};

TChain *T = 0;
  
//===================================================== Histogram Declarations
// number of histo bins
const int NTotBins = 200;
const double TotBinLow = 1.;
const double TotBinHigh = 51.;
const int RefNTotBins = 800;
const double RefTotBinLow = 1.;
const double RefTotBinHigh = 201.;

//const int num_bad = 0;
//

const int num_bad = 5;

const int bad_channels[] = {
	1161, 1472, 1670, 1696, 1896  
};

//const int num_bad = 24;
//
//const int bad_channels[] = {
//	61, 
//	64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
//	408,417,1286,
//	2124, 2404,2406,2414
//};

//const int num_bad = 63;

//const int bad_channels[] = {
//	35,40,42,45, 
//	80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,
//	136,169,183,184,192,193,194,198,199,
//	314,
//	390,391,392,395,397,433,506,507,
//	678,800,902,933,949,
//	1177,1205,1215,1255,
//	1296,1297,1299,1302,1303,1304,1307,1311,
//	1912,2140,
//	2420,2422,2430,
//	2629,2630,2662
//};
	

/*
const int bad_channels[] = {
		45,   224, 226, 234, 241, 284, 240, 604, 690, 698, 702,
		864, 1174,1437,1488,1726,1912,2131,2132,2213,2214,
		2411,2417,2425,2426,2427,2434,2436,2437,2438,2448, 2457,
		2466,2467,2468,2469,2470,2471,2472,2473,2474,2475,
		2476,2477,2478,2479,2502,2517,2518,2519,2520,2522,
		2592,2607,2613,2614,2616,2620,2624,2628,2632,2633,
		2634,2656,2661,2664,2665,2666,2667,2673,2675,2681,
		2683,2687};
*/

// Raw hits ie all hits
TH1F *hRawLe[nTdc];
TH1F *hRawTe[nTdc];
TH1F *hRawTot[nTdc];
TH1F *hGoodLe[nTdc];
TH1F *hGoodTe[nTdc];
TH1F *hGoodTot[nTdc];

TH1F *hAllRawLe;
TH1F *hAllRawTe;
TH1F *hAllRawTot;
TH1F *hAllRawPMT;
TH1F *hAllRawBar;
TH1F *hAllGoodLe;
TH1F *hAllGoodTe;
TH1F *hAllGoodTot;
TH1F *hAllGoodPMT;
TH1F *hAllGoodBar;

TH2F *h2AllGoodLe;
TH2F *h2AllGoodTe;
TH2F *h2AllGoodTot;

TH2F *h2TDCTOTvsLE;
TH2F *h2CDetX1vsX2;

TH2F *h2TOTvsXDiff1;
TH2F *h2TOTvsXDiff2;
TH2F *h2LEvsXDiff1;
TH2F *h2LEvsXDiff2;

TH2F *hBarRateHV;

TH1F *hRefRawLe;
TH1F *hRefRawTe;
TH1F *hRefRawTot;
TH1F *hRefRawPMT;
TH1F *hRefGoodLe;
TH1F *hRefGoodTe;
TH1F *hRefGoodTot;
TH1F *hRefGoodPMT;

TH1F *hMultiplicityL[nTdc];
TH1F *hMultiplicity;


// hit channel id
TH1F *hHitPMT;
TH1F *hRow;
TH1F *hRowLayer1Side1;
TH1F *hRowLayer1Side2;
TH1F *hRowLayer2Side1;
TH1F *hRowLayer2Side2;
TH1F *hLayer;
TH1F *hCol;

TH1F *hnhits1;
TH1F *hngoodhits1;
TH1F *hngoodTDChits1;
TH1F *hnhits2;
TH1F *hngoodhits2;
TH1F *hngoodTDChits2;
TH1F *hnhits_ev;
TH1F *hngoodhits_ev;
TH1F *hngoodTDChits_ev;

TH1F *hnpaddles;
TH1F *hngoodpaddles;
TH1F *hngoodTDCpaddles;

TH1F *hHitX;
TH1F *hHitY;
TH1F *hHitZ;

TH2F *hHitXY1;
TH2F *hHitXY2;

TH1F *hXECal;
TH1F *hYECal;
TH1F *hEECal;

TH2F *hXECalCDet1;
TH2F *hXECalCDet2;
TH2F *hYECalCDet1;
TH2F *hYECalCDet2;
TH2F *hEECalCDet1;
TH2F *hEECalCDet2;

TH2F *hXYECal;

TH1F *hXDiffECalCDet1;
TH1F *hXPlusECalCDet1;
TH1F *hXDiffECalCDet2;
TH1F *hXPlusECalCDet2;

TH2F *hXCDet1CDet2;
  
// 2D histograms
TH2F* h2d_RawLE;
TH2F* h2d_RawTE;
TH2F* h2d_RawTot;

TH2F* h2d_GoodLE;
TH2F* h2d_GoodTE;
TH2F* h2d_GoodTot;

TH2F* h2d_Mult;

using namespace std;

bool check_bad(int pmt, bool suppress_bad) {
	bool flag = false;
	if (!suppress_bad) return flag;
	for (int i=0;i<num_bad;i++) {
		if (pmt == bad_channels[i])  flag = true;
	}
	return flag;
}

std::vector<double> extractBinContents(const TH1* hist) {
    if (!hist) {
        throw std::invalid_argument("Null histogram pointer passed.");
    }

    int nBins = hist->GetNbinsX();
    std::vector<double> contents;

    // Loop over all *visible* bins (skip underflow bin 0 and overflow bin nBins+1)
    for (int i = 1; i <= nBins; ++i) {
        contents.push_back(hist->GetBinContent(i));
    }

    return contents;
}


vector<vector<double>> readDataFromFiles(const vector<string>& filenames) {
    const int NUM_VALUES = 42;
    vector<vector<double>> allData;

    for (const auto& filename : filenames) {
        ifstream infile(filename);
        if (!infile) {
            cerr << "Error opening file: " << filename << endl;
            continue; // Skip this file and move on
        }

        vector<double> fileData;
        double value;
        for (int i = 0; i < NUM_VALUES; ++i) {
            infile >> value;
            if (!infile) {
                cerr << "Error reading value " << i << " from file: " << filename << endl;
                break; // Stop reading this file
            }
            fileData.push_back(value);
        }

        if (fileData.size() == NUM_VALUES) {
            allData.push_back(fileData);
        } else {
            cerr << "Incomplete data in file: " << filename << endl;
        }

        infile.close();
    }

    return allData;
}

int getPixelID(int layerNum, int sideNum, int submoduleNum, int pmtNum, int pixelNum){
  // Calculate paddle number, note that missing pixels are included here
  // Validate inputs
  if (layerNum < 1 || layerNum > 2 ||
    sideNum < 1 || sideNum > 2 ||
    submoduleNum < 1 || submoduleNum > 3 ||
    pmtNum < 1 || pmtNum > 14 ||
    pixelNum < 1 || pixelNum > 16) {
    std::cerr << "Error: Invalid input values.\n"
              << "  layerNum must be 1 or 2\n"
              << "  sideNum must be 0 or 1\n"
              << "  submoduleNum must be 1 to 3\n"
              << "  pmtNum must be 1 to 14\n"
              << "  pixelNum must be 1 to 16\n";
    return -1;  // Error code
}

  int pixel = (layerNum - 1) * 1344 + //1344 pixels per layer, 0-1343 in layer 1, 1344-2687 in layer 2
  (submoduleNum - 1) * 224 + //224 pixels per side of a module
  (sideNum - 1) * 672 + //672 pixels per side
  (pmtNum-1) * 16 + //16 pixels per pmt
  pixelNum;
  return pixel - 1;
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

void PlotHVScanHighCurrent(Int_t RunNumber1=3867, Int_t nevents=50000, Int_t neventsr=500000,
	Double_t LeMin = 10.0, Double_t LeMax = 35.0,
	Double_t TotMin = 18.0, Double_t TotMax = 45.0, 
	Int_t nhitcutlow1 = 1, Int_t nhitcuthigh1 = 100,
	Int_t nhitcutlow2 = 1, Int_t nhitcuthigh2 = 100,
	Double_t XDiffCut = 0.08, Double_t XOffset = 0.02,
        Int_t layer_choice=3,	
	bool suppress_bad = false,
	Int_t nruns=30
	){

	Double_t RefLeMin = 1.0;
	Double_t RefLeMax = 251.0;
	int RefNTDCBins = (RefLeMax-RefLeMin)/4;
	Double_t RefTotMin = 1.0;
	Double_t RefTotMax = 251.0;

	int NTDCBins = 2*(LeMax-LeMin)/.0160167; // 4 ns is the trigger time, 0.018 ns is the expected time resolution, if we use a reference TDC ? 
					// 4 ns resolution is the best we can hope for, I think, using only the module trigger time.

	int NXDiffBins = (int)((2*XDiffCut)/0.0073);
	double XDiffLow = XOffset-XDiffCut;
	double XDiffHigh = XOffset+XDiffCut;

  // InFile is the input file without absolute path and without .root suffix
  // nevents is how many events to analyse, -1 for all
  
  // To execute
  // root -l
  // .L PlotRawTDC2D.C+
  // PlotRawTDC2D("filename", -1)
  double TDCBinLow = LeMin;
  double TDCBinHigh = LeMax;
  double RefTDCBinLow = RefLeMin;
  double RefTDCBinHigh = RefLeMax;

  
  // hit channel id
  
  hHitPMT = new TH1F("hHitPMT","hHitPMT",nTdc,0,nTdc);

  hnhits1 = new TH1F("hnhits1","hnhits1",150,1,151);
  hngoodhits1 = new TH1F("hngoothits1","hngoodhits1",75,1,76);
  hngoodTDChits1 = new TH1F("hngoodTDChits1","hngoodTDChits1",75,1,76);

  hnhits2 = new TH1F("hnhits2","hnhits2",150,1,151);
  hngoodhits2 = new TH1F("hngoothits2","hngoodhits2",75,1,76);
  hngoodTDChits2 = new TH1F("hngoodTDChits2","hngoodTDChits2",75,1,76);

  hnhits_ev = new TH1F("hnhits_ev","hnhits_ev",500,0,50000);
  hngoodhits_ev = new TH1F("hngoothits_ev","hngoodhits_ev",500,0,50000);
  hngoodTDChits_ev = new TH1F("hngoodTDChits_ev","hngoodTDChits_ev",500,0,50000);

  hnpaddles = new TH1F("hnpaddles","hnpaddles",200,1,201);
  hngoodpaddles = new TH1F("hngoodpaddles","hnTDCpaddles",75,1,71);
  hngoodTDCpaddles = new TH1F("hngoodTDCpaddles","hngoodTDCpaddles",75,1,76);

  hRow = new TH1F("RowNumber","RowNumber",680, 0, 680);
  hRowLayer1Side1 = new TH1F("RowNumberL1S1","RowNumberL1S1",680, 0, 680);
  hRowLayer1Side2 = new TH1F("RowNumberL1S2","RowNumberL1S2",680, 0, 680);
  hRowLayer2Side1 = new TH1F("RowNumberL2S1","RowNumberL2S1",680, 0, 680);
  hRowLayer2Side2 = new TH1F("RowNumberL2S2","RowNumberL2S2",680, 0, 680);
  hLayer = new TH1F("LayerNumber","LayerNumber",3, 0, 3);
  hCol = new TH1F("ColNumber","ColNumber",3, 0, 3);

  hHitX = new TH1F("HitXposition","HitXPosition",1000,-2.0,2.0);
  hHitY = new TH1F("HitYposition","HitYPosition",200,-0.5,0.5);
  hHitZ = new TH1F("HitZposition","HitZPosition",200,7.5,8.0);
  
  hHitXY1 = new TH2F("HitXY1position","HitXY1Position",9,-1.0,1.0,800,-2.0,2.0);
  hHitXY2 = new TH2F("HitXY2position","HitXY2Position",9,-1.0,1.0,800,-2.0,2.0);
  
  hXECal = new TH1F("XEcal","XEcal",200,-1.5,1.5);
  hYECal = new TH1F("YEcal","YEcal",200,-1.0,1.0);
  hEECal = new TH1F("YEcal","YEcal",200,0.0,20.0);
  
  hXECalCDet1 = new TH2F("XECalCDet1","XECalCDet1",100,-2.0,2.0,100,-2.0,2.0);
  hXECalCDet2 = new TH2F("XECalCDet2","XECalCDet2",100,-2.0,2.0,100,-2.0,2.0);
  hYECalCDet1 = new TH2F("YECalCDet1","YECalCDet1",100,-1.0,1.0,9,-1.0,1.0);
  hYECalCDet2 = new TH2F("YECalCDet2","YECalCDet2",100,-1.0,1.0,9,-1.0,1.0);
  hEECalCDet1 = new TH2F("EECalCDet1","EECalCDet1",100,0.0,20.0,100,-2.0,2.0);
  hEECalCDet2 = new TH2F("EECalCDet1","EECalCDet1",100,0.0,20.0,100,-2.0,2.0);
  
  hXYECal = new TH2F("XYECal","XYECal",200,-2.0,2.0,200,-2.0,2.0);
  
  hXDiffECalCDet1 = new TH1F("XDiffECalCDet1","XDiffECalCDet1",NXDiffBins,XDiffLow,XDiffHigh);
  hXPlusECalCDet1 = new TH1F("XPlusECalCDet1","XPlusECalCDet1",NXDiffBins,XDiffLow,XDiffHigh);
  hXDiffECalCDet2 = new TH1F("XDiffECalCDet2","XDiffECalCDet2",NXDiffBins,XDiffLow,XDiffHigh);
  hXPlusECalCDet2 = new TH1F("XPlusECalCDet2","XPlusECalCDet2",NXDiffBins,XDiffLow,XDiffHigh);
  
  hXCDet1CDet2 = new TH2F("XCDet1CDet2","XCDet1CDet2",200,-0.5,0.5,200,-0.5,0.5);
  
  // 2D histograms
  h2d_RawLE  = new TH2F("h2d_RawLE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTE  = new TH2F("h2d_RawTE","", NTDCBins,TDCBinLow,TDCBinHigh,nTdc+1,0,nTdc+1);
  h2d_RawTot = new TH2F("h2d_RawTot","", NTotBins,TotBinLow,TotBinHigh,nTdc+1,0,nTdc+1);
  h2d_Mult   = new TH2F("h2d_Mult","", 100,0,100,nTdc+1,0,nTdc+1);
  
  hMultiplicity = new TH1F("hMultiplicity","hMultiplicity",20,0,20);
  
  for(Int_t tdc=0; tdc<nTdc; tdc++){
    hMultiplicityL[tdc] =  new TH1F(TString::Format("hMultiplicity_Bar%d",tdc),
		      TString::Format("hMultiplicity_Bar%d",tdc),
		      10, 0, 10);
  }// element loop

  hAllRawLe = new TH1F(TString::Format("hRawLe"),
            TString::Format("hRawLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllRawTe = new TH1F(TString::Format("hRawTe"),
            TString::Format("hRawTe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllRawTot = new TH1F(TString::Format("hRawTot"),
            TString::Format("hRawTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllRawPMT = new TH1F(TString::Format("hRawPMT"),
            TString::Format("hRawPMT"),
            nTdc, 0, nTdc);
  hAllRawBar = new TH1F(TString::Format("hRawBar"),
            TString::Format("hRawBar"),
            168, 0, 168);
  hAllGoodLe = new TH1F(TString::Format("hAllGoodLe"),
            TString::Format("hAllGoodLe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllGoodTe = new TH1F(TString::Format("hAllGoodTe"),
            TString::Format("hAllGoodTe"),
            NTDCBins, TDCBinLow, TDCBinHigh);
  hAllGoodTot = new TH1F(TString::Format("hAllGoodTot"),
            TString::Format("hAllGoodTot"),
            NTotBins, TotBinLow, TotBinHigh);
  hAllGoodPMT = new TH1F(TString::Format("hAllGoodPMT"),
            TString::Format("hAllGoodPMT"),
            nTdc, 0, nTdc);
  hAllGoodBar = new TH1F(TString::Format("hAllGoodBar"),
            TString::Format("hAllGoodBar"),
            168, 0, 168);
  h2AllGoodLe = new TH2F(TString::Format("h2AllGoodLe"),
            TString::Format("h2AllGoodLe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTe = new TH2F(TString::Format("h2AllGoodTe"),
            TString::Format("h2AllGoodTe"),nTdc,0,nTdc,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2AllGoodTot = new TH2F(TString::Format("h2AllGoodTot"),
            TString::Format("h2AllGoodTot"),nTdc,0,nTdc,
            NTotBins, TotBinLow, TotBinHigh);
  hBarRateHV = new TH2F(TString::Format("hBarRateHV"),
	    TString::Format("hBarRateHV"), 250,1,250000,
	    50,600,800);

  h2TDCTOTvsLE = new TH2F(TString::Format("h2TDCTOTvsLE"),
            TString::Format("h2TDCTOTvsLE"),NTotBins,TotBinLow,TotBinHigh,
            NTDCBins, TDCBinLow, TDCBinHigh);
  h2CDetX1vsX2 = new TH2F(TString::Format("h2CDetX1vsX2"),
            TString::Format("h2CDetX1vsX2"),1000, -2.0, 2.0, 
            1000, -2.0, 2.0);
  h2TOTvsXDiff1 = new TH2F(TString::Format("h2TOTvsXDiff1"),
            TString::Format("h2TOTvsXDiff1"),NTotBins,TotBinLow,TotBinHigh,
            1000, -0.3, 0.3);
  h2TOTvsXDiff2 = new TH2F(TString::Format("h2TOTvsXDiff2"),
            TString::Format("h2TOTvsXDiff2"),NTotBins,TotBinLow,TotBinHigh,
            1000, -0.3, 0.3);
  h2LEvsXDiff1 = new TH2F(TString::Format("h2LEvsXDiff1"),
            TString::Format("h2LEvsXDiff1"),NTDCBins,TDCBinLow,TDCBinHigh,
            1000, -0.3, 0.3);
  h2LEvsXDiff2 = new TH2F(TString::Format("h2LEvsXDiff2"),
            TString::Format("h2LEvsXDiff2"),NTDCBins,TDCBinLow,TDCBinHigh,
            1000, -0.3, 0.3);

  hRefRawLe = new TH1F(TString::Format("hRefRawLe"),
            TString::Format("hRefRawLe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefRawTe = new TH1F(TString::Format("hRefRawTe"),
            TString::Format("hRefRawTe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefRawTot = new TH1F(TString::Format("hRefRawTot"),
            TString::Format("hRefRawTot"),
            RefNTotBins, RefTotBinLow, RefTotBinHigh);
  hRefRawPMT = new TH1F(TString::Format("hRefRawPMT"),
            TString::Format("hRefRawPMT"),
            32, 2688, 2720);
  hRefGoodLe = new TH1F(TString::Format("hRefGoodLe"),
            TString::Format("hRefGoodLe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTe = new TH1F(TString::Format("hRefGoodTe"),
            TString::Format("hRefGoodTe"),
            RefNTDCBins, RefTDCBinLow, RefTDCBinHigh);
  hRefGoodTot = new TH1F(TString::Format("hRefGoodTot"),
            TString::Format("hRefGoodTot"),
            RefNTotBins, RefTotBinLow, RefTotBinHigh);
  hRefGoodPMT = new TH1F(TString::Format("hRefGoodPMT"),
            TString::Format("hRefGoodPMT"),
            2720, 0, 2720);
  
  
  for(Int_t bar=0; bar<(nTdc); bar++){
    // raw hits
    // leading edge
    hRawLe[bar] = new TH1F(TString::Format("hRawLe_Bar%d",bar),
 	    TString::Format("hRawLe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // trailing edge 
    hRawTe[bar] = new TH1F(TString::Format("hRawTe_Bar%d",bar),
	    TString::Format("hRawTe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // tot 
    hRawTot[bar] = new TH1F(TString::Format("hRawTot_Bar%d",bar),
			     TString::Format("hRawTot_Bar%d",bar),
			     NTotBins, TotBinLow, TotBinHigh);
    // good hits
    // leading edge
    hGoodLe[bar] = new TH1F(TString::Format("hGoodLe_Bar%d",bar),
 	    TString::Format("hGoodLe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // trailing edge 
    hGoodTe[bar] = new TH1F(TString::Format("hGoodTe_Bar%d",bar),
	    TString::Format("hGoodTe_Bar%d",bar),
	    NTDCBins, TDCBinLow, TDCBinHigh);
    // tot 
    hGoodTot[bar] = new TH1F(TString::Format("hGoodTot_Bar%d",bar),
			     TString::Format("hGoodTot_Bar%d",bar),
			     NTotBins, TotBinLow, TotBinHigh);
  }// bar loop



  //========================================================= Get data from tree
  if(!T) { 
    // TString sInFile = REPLAYED_DIR + "/" + InFile + ".root";
    T = new TChain("T");

    TString subfile, sInFile;

    subfile = TString::Format("cdet_%d_%d",RunNumber1,neventsr);
    sInFile = REPLAYED_DIR + "/" + subfile + ".root";
    cout << "Input ROOT file = " << sInFile << endl;
    cout << "Adding " << sInFile << endl;
    T->Add(sInFile);
    cout << "Adding " << nruns << " files ... " << endl;
    for (Int_t i=1; i<=nruns; i++) {
        subfile = TString::Format("cdet_%d_%d_%d",RunNumber1,neventsr,i);
        //subfile = TString::Format("_%d_1000000_%d",RunNumber1,i);
        sInFile = REPLAYED_DIR + "/" + subfile + ".root";
        cout << "Input ROOT file = " << sInFile << endl;
        cout << "Adding " << sInFile << endl;
        T->Add(sInFile);
    }
 
    
    // disable all branches
    T->SetBranchStatus("*",0);
    // enable branches
    T->SetBranchStatus("earm.cdet.*",1);
    T->SetBranchStatus("earm.ecal.*",1);

    T->SetBranchAddress("earm.cdet.tdc_mult",TCDet::TDCmult);
    
    T->SetBranchAddress("earm.cdet.hits.TDCelemID",TCDet::RawElID);
    T->SetBranchAddress("earm.cdet.hits.t",TCDet::RawElLE);
    T->SetBranchAddress("earm.cdet.hits.t_te",TCDet::RawElTE);
    T->SetBranchAddress("earm.cdet.hits.t_tot",TCDet::RawElTot);

    T->SetBranchAddress("earm.cdet.hit.pmtnum",TCDet::GoodElID);
    T->SetBranchAddress("earm.cdet.hit.tdc_le",TCDet::GoodElLE);
    T->SetBranchAddress("earm.cdet.hit.tdc_te",TCDet::GoodElTE);
    T->SetBranchAddress("earm.cdet.hit.tdc_tot",TCDet::GoodElTot);
    
    T->SetBranchAddress("earm.cdet.hit.xhit",TCDet::GoodX);
    T->SetBranchAddress("earm.cdet.hit.yhit",TCDet::GoodY);
    T->SetBranchAddress("earm.cdet.hit.zhit",TCDet::GoodZ);
    
    T->SetBranchAddress("earm.cdet.hit.row",TCDet::GoodCol);
    T->SetBranchAddress("earm.cdet.hit.col",TCDet::GoodRow);
    T->SetBranchAddress("earm.cdet.hit.layer",TCDet::GoodLayer);

        
    T->SetBranchAddress("earm.cdet.nhits",&TCDet::nhits);
    T->SetBranchAddress("earm.cdet.ngoodhits",&TCDet::ngoodhits);
    T->SetBranchAddress("earm.cdet.ngoodTDChits",&TCDet::ngoodTDChits);
    T->SetBranchAddress("earm.ecal.x",&TCDet::GoodECalX);
    T->SetBranchAddress("earm.ecal.y",&TCDet::GoodECalY);
    T->SetBranchAddress("earm.ecal.e",&TCDet::GoodECalE);

    // enable vector size branches
    T->SetBranchAddress("Ndata.earm.cdet.tdc_mult",&TCDet::NdataMult); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.TDCelemID",&TCDet::NdataRawElID); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.t",&TCDet::NdataRawElLE); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.t_te",&TCDet::NdataRawElTE); 
    T->SetBranchAddress("Ndata.earm.cdet.hits.t_tot",&TCDet::NdataRawElTot); 
    
    T->SetBranchAddress("Ndata.earm.cdet.hit.pmtnum",&TCDet::NdataGoodElID);
    T->SetBranchAddress("Ndata.earm.cdet.hit.tdc_le",&TCDet::NdataGoodElLE);
    T->SetBranchAddress("Ndata.earm.cdet.hit.tdc_te",&TCDet::NdataGoodElTE);
    T->SetBranchAddress("Ndata.earm.cdet.hit.tdc_tot",&TCDet::NdataGoodElTot);
    
    T->SetBranchAddress("Ndata.earm.cdet.hit.xhit",&TCDet::NdataGoodX);
    T->SetBranchAddress("Ndata.earm.cdet.hit.yhit",&TCDet::NdataGoodY);
    T->SetBranchAddress("Ndata.earm.cdet.hit.zhit",&TCDet::NdataGoodZ);
    
    T->SetBranchAddress("Ndata.earm.cdet.hit.row",&TCDet::NdataGoodCol);
    T->SetBranchAddress("Ndata.earm.cdet.hit.col",&TCDet::NdataGoodRow);
    T->SetBranchAddress("Ndata.earm.cdet.hit.layer",&TCDet::NdataGoodLayer);
    

  }//setting tree
  
  //========================================================= Check no of events
  Int_t Nev = T->GetEntries();
  cout << "N entries in tree is " << Nev << endl;
  Int_t NEventsAnalysis;
  if(nevents==-1) NEventsAnalysis = Nev;
  else NEventsAnalysis = nevents;
  cout << "Running analysis for " << NEventsAnalysis << " events" << endl;
  

  
  //==================================================== Create output root file
  // root file for viewing fits
  TString subfile; 
  subfile = TString::Format("gep5_replayed_nogems_%d_50k_events.root",RunNumber1);
  TString outrootfile = ANALYSED_DIR + "/RawTDC_" + subfile;
  TFile *f = new TFile(outrootfile, "RECREATE");



  //================================================================= Event Loop
  // variables outside event loop
  Int_t EventCounter = 0;
  cout << "Starting Event Loop" << endl;

    int eff_denominator = 0;
    int eff_numerator_layer1 = 0;
    int eff_numerator_layer2 = 0;
    int eff_numerator = 0;

  // event loop start
  for(Int_t event=0; event<NEventsAnalysis; event++){
    
    T->GetEntry(event);
    EventCounter++;
    if (EventCounter % 1000 == 0) {
	cout << EventCounter << "/" << NEventsAnalysis << "/ Nhits = " << (Int_t)TCDet::nhits << endl;
    	for (Int_t nfill=0; nfill<(Int_t)TCDet::nhits; nfill++) {hnhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<(Int_t)TCDet::ngoodhits; nfill++) {hngoodhits_ev->Fill(EventCounter);}
    	for (Int_t nfill=0; nfill<(Int_t)TCDet::ngoodTDChits; nfill++) {hngoodTDChits_ev->Fill(EventCounter);}
    }


    //cout << "Raw TDC hit loop: " << TCDet::NdataRawElID << endl;
    
    // First pass through hits:  purpose is to get reference LE TDC Value for this event
    
    double event_ref_tdc = 0.0;

    for(Int_t el=0; el<TCDet::NdataRawElID; el++) {
 
	bool good_ref_le_time = TCDet::RawElLE[el] >= 0.0/TDC_calib_to_ns && TCDet::RawElLE[el] <= 100.0/TDC_calib_to_ns;
	bool good_ref_tot = TCDet::GoodElTot[el] >= 0.0/TDC_calib_to_ns && TCDet::GoodElTot[el] <= 200.0/TDC_calib_to_ns;
	bool good_ref_event = good_ref_le_time && good_ref_tot;
	if ( good_ref_event ) {

	  //if (TCDet::RawElID[el] > 2687) {
	  //	cout << "el = " << el << " Raw ID = " << TCDet::RawElID[el] << " raw le = " << 
	//	TCDet::RawElLE[el] << " raw te = " << TCDet::RawElTE[el] << " raw tot = " << 
	//	TCDet::RawElTot[el] << " CDet X = " << TCDet::GoodX[el] << " ECal X = " << TCDet::GoodECalX << endl;
	  //}
	  if ( !check_bad(TCDet::RawElID[el],suppress_bad) ) {
	   //cout << " el = " << el << endl;
	   //cout << " tdc = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
	   if ( (Int_t)TCDet::RawElID[el] == 2696 && (Int_t)TCDet::RawElLE[el]>0 && (Int_t)TCDet::RawElTot[el]>0 ) {
	    //cout << " Ref  ID = " << (Int_t)TCDet::RawElID[el] << " el = " << el << "    LE = " << TCDet::RawElLE[el]*TDC_calib_to_ns 
	    //		<< "    TE = " << TCDet::RawElTE[el]*TDC_calib_to_ns << "    ToT = " << TCDet::RawElTot[el]*TDC_calib_to_ns << endl;
 	    hRefRawLe->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns);
	    hRefRawTe->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns);
	    hRefRawTot->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns);
	    hRefRawPMT->Fill(TCDet::RawElID[el]);

	    event_ref_tdc = TCDet::RawElLE[el]*TDC_calib_to_ns;
	

	   }
	  }
	}
    }// end ref TDC loop
    

    // second pass: fill raw CDet TDC histos
    
    for(Int_t el=0; el<TCDet::NdataRawElID; el++){

	bool good_raw_le_time = TCDet::RawElLE[el] >= LeMin/TDC_calib_to_ns && TCDet::RawElLE[el] <= LeMax/TDC_calib_to_ns;
	bool good_raw_tot = TCDet::GoodElTot[el] >= TotMin/TDC_calib_to_ns && TCDet::GoodElTot[el] <= TotMax/TDC_calib_to_ns;
	bool good_mult = TCDet::TDCmult[el] < TDCmult_cut;


	bool good_raw_event = good_raw_le_time && good_raw_tot && good_mult;


	//if ((Int_t)TCDet::RawElID[el] > 1000) cout << "el = " << el << " Hit ID = " << (Int_t)TCDet::RawElID[el] << "    TDC = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
	//cout << "Raw ID = " << TCDet::RawElID[el] << " raw le = " << TCDet::RawElLE[el] << " raw te = " << TCDet::RawElTE[el] << " raw tot = " << TCDet::RawElTot[el] << endl;
	if ( good_raw_event ) {

	  //if (TCDet::RawElID[el] > 2687) {
	  //	cout << "el = " << el << " Raw ID = " << TCDet::RawElID[el] << " raw le = " << 
	//	TCDet::RawElLE[el] << " raw te = " << TCDet::RawElTE[el] << " raw tot = " << 
	//	TCDet::RawElTot[el] << " CDet X = " << TCDet::GoodX[el] << " ECal X = " << TCDet::GoodECalX << endl;
	  //}
	  if ( !check_bad(TCDet::RawElID[el],suppress_bad) ) {
	   //cout << " el = " << el << endl;
	   //cout << " tdc = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
	   if ( (Int_t)TCDet::RawElID[el] < 2688 ) {
	    //if ((Int_t)TCDet::RawElID[el] > nTdc) cout << " CDet ID = " << (Int_t)TCDet::RawElID[el] << "    TDC = " << TCDet::RawElLE[el]*TDC_calib_to_ns << endl;
	    hRawLe[(Int_t)TCDet::RawElID[el]]->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    hRawTe[(Int_t)TCDet::RawElID[el]]->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    hRawTot[(Int_t)TCDet::RawElID[el]]->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns);
 	    hAllRawLe->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    hAllRawTe->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    hAllRawTot->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns);
	    hAllRawPMT->Fill(TCDet::RawElID[el]);
	    hAllRawBar->Fill((Int_t)(TCDet::RawElID[el]/16));

	    h2d_RawLE->Fill(TCDet::RawElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0, (Int_t)TCDet::RawElID[el]);
	    h2d_RawTE->Fill(TCDet::RawElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0, (Int_t)TCDet::RawElID[el]);
	    h2d_RawTot->Fill(TCDet::RawElTot[el]*TDC_calib_to_ns, (Int_t)TCDet::RawElID[el]);

	   }
	  }
	}
	

    }// all raw tdc hit loop


    // Third pass:  Get layer occupancies
    
    int nhitsc1 = 0;
    int nhitsc2 = 0;
    int ngoodhitsc1 = 0;
    int ngoodhitsc2 = 0;
    int ngoodTDChitsc1 = 0;
    int ngoodTDChitsc2 = 0;
    for (int j=0; j<nTdc; j++) {
	TCDet::nhits_paddles[j]=0;
	TCDet::ngoodhits_paddles[j]=0;
	TCDet::ngoodTDChits_paddles[j]=0;
    }
    TCDet::npaddles=0;
    TCDet::ngoodpaddles=0;
    TCDet::ngoodTDCpaddles=0;
    
    for(Int_t el=0; el<TCDet::NdataGoodElID; el++){
	int sbselemid = (Int_t)TCDet::GoodElID[el];
	int sbsrown = sbselemid%672;
	int sbscoln = sbselemid/672;
	//int sbsrown = (Int_t)TCDet::GoodRow[el];
	//int sbscoln = (Int_t)TCDet::GoodCol[el];
	int mylayern = sbscoln/2;
	int mypaddlen = sbscoln*672 + sbsrown;

	if (mylayern == 0) {
		nhitsc1++;
	} else {
		nhitsc2++;
	}
	TCDet::nhits_paddles[mypaddlen]++;

        bool good_ecal_reconstruction = TCDet::GoodECalY > -1.2 && TCDet::GoodECalY < 1.2 &&
			TCDet::GoodECalX > -1.5 && TCDet::GoodECalX < 1.5 &&
			TCDet::GoodECalX != 0.00 && TCDet::GoodECalY != 0.00 ;
	bool good_le_time = TCDet::GoodElLE[el] >= LeMin/TDC_calib_to_ns && TCDet::GoodElLE[el] <= LeMax/TDC_calib_to_ns;
	bool good_tot = TCDet::GoodElTot[el] >= TotMin/TDC_calib_to_ns && TCDet::GoodElTot[el] <= TotMax/TDC_calib_to_ns;
	bool good_hit_mult = TCDet::TDCmult[el] < TDCmult_cut;
	bool good_cdet_X = TCDet::GoodX[el] < xcut;
	bool good_ecal_diff_x = (TCDet::GoodX[el]-(TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)-XOffset) <= XDiffCut && 
			(TCDet::GoodX[el]-(TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)-XOffset) >= -1.0*XDiffCut;
	bool good_ecal_diff_y = (TCDet::GoodY[el]-(TCDet::GoodECalY*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)) <= cdet_y_half_length && 
			(TCDet::GoodY[el]-(TCDet::GoodECalY*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)) >= -1.0*cdet_y_half_length;
	

	bool good_CDet_event = good_ecal_reconstruction && good_ecal_diff_x && good_ecal_diff_y && good_le_time && good_tot && good_hit_mult && good_cdet_X;


	if (good_CDet_event) {

	  if ( !check_bad(TCDet::GoodElID[el], suppress_bad) ) {
	   if ( (Int_t)TCDet::GoodElID[el]%NumSidesTotal < NumCDetPaddlesPerSide )  {

	    //cout << "Hit number " << el << ":    Paddle = " << mypaddlen << " Row = " << sbsrown  << " Col = " << sbscoln  << " hits = " << TCDet::ngoodTDChits_paddles[mypaddlen] << endl;
	    //cout << "el = " << el << " Good ID = " << TCDet::GoodElID[el] << " Good le = " << 
	//	TCDet::GoodElLE[el] << " Good te = " << TCDet::GoodElTE[el] << " Good tot = " << 
	//	TCDet::GoodElTot[el] << " CDet X = " << TCDet::GoodX[el] << " ECal X = " << TCDet::GoodECalX << endl;
	    if (mylayern == 0) {
	        ngoodhitsc1++;
	    } else {
	        ngoodhitsc2++;
	    }
	    TCDet::ngoodhits_paddles[mypaddlen]++;
	   }
	  }
	}
    }
    for (int j=0; j<nTdc; j++) {
	if (TCDet::nhits_paddles[j] > 0) {
		 TCDet::npaddles++;
		 //cout << "Paddle = " << j <<  "  nhits = " << TCDet::ngoodTDChits_paddles[j] << endl;
	}
    }
    for (int j=0; j<nTdc; j++) {
	if (TCDet::ngoodhits_paddles[j] > 0) {
		 TCDet::ngoodpaddles++;
		 //cout << "Paddle = " << j <<  "  nhits = " << TCDet::ngoodTDChits_paddles[j] << endl;
	}
    }
    //cout << "event " << event << endl;
    //cout << "Number of good layer 1 hits: " << ngoodTDChitsc1 << endl;
    //cout << "Number of good layer 2 hits: " << ngoodTDChitsc2 << endl;
    //cout << "Layer 1 Hit Cut " << nhitcutlow1 << " " << nhitcuthigh1 << endl;
    //cout << "Layer 2 Hit Cut " << nhitcutlow2 << " " << nhitcuthigh2 << endl;

    hnpaddles->Fill(TCDet::npaddles);
    hngoodpaddles->Fill(TCDet::ngoodpaddles);

    // Fourth pass:  use layer occupancies to apply additional cuts
 
    for(Int_t el=0; el<TCDet::NdataGoodElID; el++){
        bool goodhit_ecal_reconstruction = TCDet::GoodECalY > -1.2 && TCDet::GoodECalY < 1.2 &&
			TCDet::GoodECalX > -1.5 && TCDet::GoodECalX < 1.5 &&
			TCDet::GoodECalX != 0.00 && TCDet::GoodECalY != 0.00 ;
	bool goodhit_le_time = TCDet::GoodElLE[el] >= LeMin/TDC_calib_to_ns && TCDet::GoodElLE[el] <= LeMax/TDC_calib_to_ns;
	bool goodhit_tot = TCDet::GoodElTot[el] >= TotMin/TDC_calib_to_ns && TCDet::GoodElTot[el] <= TotMax/TDC_calib_to_ns;
	bool goodhit_hit_mult = TCDet::TDCmult[el] < TDCmult_cut;
	bool goodhit_cdet_X = TCDet::GoodX[el] < xcut;
	bool goodhit_low = ngoodhitsc1 >= nhitcutlow1  && ngoodhitsc2 >= nhitcutlow2;
	bool goodhit_high  = ngoodhitsc1 <= nhitcuthigh1 && ngoodhitsc2 <= nhitcuthigh2; 
	bool goodhit_ecal_diff_x = (TCDet::GoodX[el]-(TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)-XOffset) <= XDiffCut && 
			(TCDet::GoodX[el]-(TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)-XOffset) >= -1.0*XDiffCut;
	bool goodhit_ecal_diff_y = (TCDet::GoodY[el]-(TCDet::GoodECalY*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)) <= cdet_y_half_length && 
			(TCDet::GoodY[el]-(TCDet::GoodECalY*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist)) >= -1.0*cdet_y_half_length;


	bool goodhit_CDet_event = goodhit_ecal_reconstruction && goodhit_ecal_diff_x && goodhit_ecal_diff_y && goodhit_le_time && goodhit_tot 
		&& goodhit_hit_mult && goodhit_cdet_X && goodhit_low && goodhit_high;

	if (goodhit_CDet_event) {

	  if ( !check_bad(TCDet::GoodElID[el], suppress_bad) ) {
	   if ( (Int_t)TCDet::GoodElID[el]%NumSidesTotal < NumCDetPaddlesPerSide )  {
    	    //cout << "event " << event << endl;
	    //cout << "el = " << el << " Good ID = " << TCDet::GoodElID[el] << " Good le = " << 
		//TCDet::GoodElLE[el] << " Good te = " << TCDet::GoodElTE[el] << " Good tot = " << 
		//TCDet::GoodElTot[el] << " CDet X = " << TCDet::GoodX[el] << " ECal X = " << TCDet::GoodECalX << endl;

	    //cout << "Filling good timing histos ... " << ngoodTDChitsc1 << " " << endl;
	    
	    //std::cout << "Layer = " << (Int_t)TCDet::GoodLayer[el] << " Side = " << (Int_t)TCDet::GoodCol[el] << std::endl;
	    
	    int sbselem = (Int_t)TCDet::GoodElID[el];
	    int sbsrow = sbselem%672;
	    int sbscol = sbselem/672;
	    //int sbscol = (Int_t)TCDet::GoodCol[el];
	    //int sbsrow = (Int_t)TCDet::GoodRow[el];
	    int myside = sbscol%2;
	    int mylayer = sbscol/2;
	    int mypaddle = sbscol*672 + sbsrow;

	    if (mylayer == 0) {
		ngoodTDChitsc1++;
		eff_numerator_layer1++;
	    } else {
		ngoodTDChitsc2++;
		eff_numerator_layer2++;
	    }

	    TCDet::ngoodTDChits_paddles[mypaddle]++;
	    
	    if ( (layer_choice == 1 && mylayer == 0) || (layer_choice == 2 && mylayer == 1) || 
			(layer_choice == 3 && ngoodhitsc1>=1 && ngoodhitsc2 >= 1) )  
		{
		eff_numerator++;
	    	hGoodLe[(Int_t)TCDet::GoodElID[el]]->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    	hGoodTe[(Int_t)TCDet::GoodElID[el]]->Fill(TCDet::GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    	hGoodTot[(Int_t)TCDet::GoodElID[el]]->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns);
	    	hAllGoodLe->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    	hAllGoodTe->Fill(TCDet::GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    	hAllGoodTot->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns);
	    	hAllGoodPMT->Fill(TCDet::GoodElID[el]);
	    	hAllGoodBar->Fill((Int_t)(TCDet::GoodElID[el]/16));
	    	h2AllGoodLe->Fill(TCDet::GoodElID[el],TCDet::GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    	h2AllGoodTe->Fill(TCDet::GoodElID[el],TCDet::GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    	h2AllGoodTot->Fill(TCDet::GoodElID[el],TCDet::GoodElTot[el]*TDC_calib_to_ns);

	    	h2TDCTOTvsLE->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns,TCDet::GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	  
	    	hHitPMT->Fill((Int_t)TCDet::GoodElID[el]);
	    	hRow->Fill((Int_t)TCDet::GoodRow[el]);
	    }
	
	    if (myside == 0) {
		if (mylayer == 0) {
			hRowLayer1Side1->Fill((Int_t)TCDet::GoodRow[el]);
		} else {
			hRowLayer2Side1->Fill((Int_t)TCDet::GoodRow[el]);
		}
	    } else {
		if(mylayer == 0) {
			hRowLayer1Side2->Fill((Int_t)TCDet::GoodRow[el]);
		} else {
			hRowLayer2Side2->Fill((Int_t)TCDet::GoodRow[el]);
		}
	    }	
			
	    hCol->Fill(myside);
	    hLayer->Fill(mylayer);

	    hHitX->Fill(TCDet::GoodX[el]);
	    hHitY->Fill(TCDet::GoodY[el]);
	    hHitZ->Fill(TCDet::GoodZ[el]);
	    if (mylayer==0) {
	    	h2TOTvsXDiff1->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns,TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    	h2LEvsXDiff1->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0,TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
		hHitXY1->Fill(TCDet::GoodY[el],TCDet::GoodX[el]);
		hXECalCDet1->Fill(TCDet::GoodX[el],TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
		hYECalCDet1->Fill(TCDet::GoodY[el],TCDet::GoodECalY*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    	hXDiffECalCDet1->Fill(TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    	hXPlusECalCDet1->Fill(TCDet::GoodX[el]+TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
		hEECalCDet1->Fill(TCDet::GoodECalE,TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    } else {
	    	h2TOTvsXDiff2->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns,TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    	h2LEvsXDiff2->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0,TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
		hHitXY2->Fill(TCDet::GoodY[el],TCDet::GoodX[el]);
		hXECalCDet2->Fill(TCDet::GoodX[el],TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
		hYECalCDet2->Fill(TCDet::GoodY[el],TCDet::GoodECalY*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    	hXDiffECalCDet2->Fill(TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    	hXPlusECalCDet2->Fill(TCDet::GoodX[el]+TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
		hEECalCDet2->Fill(TCDet::GoodECalE,TCDet::GoodX[el]-TCDet::GoodECalX*(TCDet::GoodZ[el]-cdet_dist_offset)/ecal_dist);
	    }


	   } else {
	    hRefGoodLe->Fill(TCDet::GoodElLE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    hRefGoodTe->Fill(TCDet::GoodElTE[el]*TDC_calib_to_ns-event_ref_tdc+60.0);
	    hRefGoodTot->Fill(TCDet::GoodElTot[el]*TDC_calib_to_ns);
	    hRefGoodPMT->Fill(TCDet::GoodElID[el]*TDC_calib_to_ns);
	   }
	  }
	}


    }// all good tdc hit loop

    
    if (TCDet::GoodECalX != 0.00 && TCDet::GoodECalY != 0.00) {
      eff_denominator++;
      hXYECal->Fill(TCDet::GoodECalY,TCDet::GoodECalX);
      hXECal->Fill(TCDet::GoodECalX);
      hYECal->Fill(TCDet::GoodECalY);
      hEECal->Fill(TCDet::GoodECalE);
    };
    
    
    hnhits1->Fill(nhitsc1);
    hngoodhits1->Fill(ngoodhitsc1);
    hngoodTDChits1->Fill(ngoodTDChitsc1);
    hnhits2->Fill(nhitsc2);
    hngoodhits2->Fill(ngoodhitsc2);
    hngoodTDChits2->Fill(ngoodTDChitsc2);
    for (int j=0; j<nTdc; j++) {
	if (TCDet::ngoodTDChits_paddles[j] > 0) {
		 TCDet::ngoodTDCpaddles++;
		 //cout << "Paddle = " << j <<  "  nhits = " << TCDet::ngoodTDChits_paddles[j] << endl;
	}
    }
    hngoodTDCpaddles->Fill(TCDet::ngoodTDCpaddles);

    //cout << "Element loop: " << TCDet::NdataMult << endl;
    for(Int_t tdc=0; tdc<TCDet::NdataMult; tdc++){
      if (!check_bad(TCDet::RawElID[tdc],suppress_bad)) {
      hMultiplicity->Fill(TCDet::TDCmult[tdc]);
	hMultiplicityL[(Int_t)TCDet::RawElID[tdc]]->Fill(TCDet::TDCmult[tdc]);
	if( TCDet::TDCmult[tdc] != 0 )
	  h2d_Mult->Fill(TCDet::TDCmult[tdc], (Int_t)TCDet::RawElID[tdc] );
      }
    }// element loop


  }// event loop

  std::cout << "Candidate Events = " << eff_denominator << std::endl;
  std::cout << "Layer 1 Events = " << eff_numerator_layer1 << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator_layer1/eff_denominator  << std::endl;
  std::cout << "Layer 2 Events = " << eff_numerator_layer2 << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator_layer2/eff_denominator <<  std::endl;
  std::cout << "One Good Layer Events = " << eff_numerator << "     Avg Hits Per Candidate Event = " << 1.0*eff_numerator/eff_denominator <<  std::endl;


  for (Int_t b=0; b<NumCDetPaddles; b++) {
	//if (hRawLe[b]->GetEntries() > EventCounter/HotChannelRatio) {
	if (hRawLe[b]->GetEntries() > EventCounter/NumCDetPaddles*2*1000) {
		int myhotlayer = b/1344 + 1;
		int myhotside = (b%1344)/672 + 1;
		int myhotmodule = (b%672)/244 + 1;
		int myhotbar = (b%672)%224/16 + 1;
		int myhotpaddle = ((b%672)%224)%16 + 1;
		int mycable = b/16;
		std::cout << "Hot PMT!! ID = " << b << "  layer = " << myhotlayer <<
		"   side = " << myhotside << "   module = " << myhotmodule <<
		"   bar = " << myhotbar << "   paddle_PMT = " << myhotpaddle << " CDet Cable = " << mycable << "   Entries = " << hRawLe[b]->GetEntries() << std::endl;
	}
  }
  
    
    //========================================================== Write histos
  for(Int_t b=0; b<nTdc; b++){
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawLe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawLe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawLe[b]->Write();
    // hRawLe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTe[b]->GetXaxis()->SetTitle("time (ns)");
    // hRawTe[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTe[b]->Write();
    // hRawTe[b]->GetXaxis()->SetLabelSize(0.06);
    hRawTot[b]->GetXaxis()->SetTitle("tot (ns)");
    // hRawTot[b]->GetXaxis()->SetTitleSize(0.05);
    hRawTot[b]->Write();
    // hMultiplicityL[b]->GetXaxis()->SetLabelSize(0.06);
    hMultiplicityL[b]->GetXaxis()->SetTitle("tdc ref hit mult");
    hMultiplicityL[b]->Write();
  }
  // hMultiplicity->GetXaxis()->SetLabelSize(0.06);
  hMultiplicity->GetXaxis()->SetTitle("tdc hit multiplicity");
  hMultiplicity->Write();
  // hHitPMT->GetXaxis()->SetLabelSize(0.06);
  hHitPMT->GetXaxis()->SetTitle("Bar ID of Left PMT Hit");
  hHitPMT->SetTitle("");
  hHitPMT->Write();

  // 2D histograms
  h2d_RawLE->GetXaxis()->SetTitle("TDC Leading Edge Time [ns]");
  h2d_RawLE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawLE->SetTitle("");
  h2d_RawLE->Write();
  h2d_RawTE->GetXaxis()->SetTitle("TDC Trailing Edge Time [ns]");
  h2d_RawTE->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTE->SetTitle("");
  h2d_RawTE->Write();

  h2d_RawTot->GetXaxis()->SetTitle("TDC Time-over-threshold [ns]");
  h2d_RawTot->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_RawTot->SetTitle("");
  h2d_RawTot->Write();

  h2d_Mult->GetXaxis()->SetTitle("TDC Multiplicity [ns]");
  h2d_Mult->GetYaxis()->SetTitle("PMT number (Left)");
  h2d_Mult->SetTitle("");
  h2d_Mult->Write();

  // Get HV values
  vector<string> HVfilenames = {"l1Left.dat", "l1Right.dat", "l2Left.dat", "l2Right.dat"};

  vector<vector<double>> data = readDataFromFiles(HVfilenames);
  std::vector<double> barRateContents = extractBinContents(hAllRawBar);
  Double_t xval,yval;
  for (int ii=0;ii<4;ii++) {
	for (int jj=0;jj<42;jj++) {
		xval = data[ii][jj];
	        yval = barRateContents[ii*42+jj];
		//cout << "Contents:  " << xval << " "  << yval << endl;
		hBarRateHV->Fill(yval,-xval);
	}
  }

  //========================================================== Close output file
  f->Close();



  //================================================================== End Macro
}// end main

TCanvas *plotBarRateHV() {
  
  TCanvas *daa = new TCanvas("All TDC", "All TDC", 50,50,800,800);

  daa->cd();
  gPad->SetLogx();
  hBarRateHV->Draw();

  return daa;
}
	

TCanvas *plotAllTDC(){

  TCanvas *caa = new TCanvas("All TDC", "All TDC", 50,50,800,800);
  caa->Divide(2,3,0.01,0.01,0);
  TCanvas *caaa = new TCanvas("All Chan", "All Chan", 850,50,800,800);
  caaa->Divide(2,2,0.01,0.01,0);

  caa->cd(1);
  hAllRawLe->SetFillColor(kRed);
  hAllRawLe->SetMinimum(0.0);
  hAllRawLe->Draw();
  
  caa->cd(2);
  hAllGoodLe->SetFillColor(kBlue);
  hAllGoodLe->SetMinimum(0.0);
  hAllGoodLe->Draw();
  
  caa->cd(3);
  hAllRawTe->SetFillColor(kRed);
  hAllRawTe->SetMinimum(0.0);
  hAllRawTe->Draw();
  
  caa->cd(4);
  hAllGoodTe->SetFillColor(kBlue);
  hAllGoodTe->SetMinimum(0.0);
  hAllGoodTe->Draw();
  

  caa->cd(5);
  //gPad->SetLogy();
  hAllRawTot->SetFillColor(kRed);
  hAllRawTot->Draw();

  caa->cd(6);
  //gPad->SetLogy();
  hAllGoodTot->SetFillColor(kBlue);
  hAllGoodTot->Draw();

  caaa->cd(1);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllRawPMT->SetFillColor(kRed);
  hAllRawPMT->Draw();
  
  caaa->cd(2);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllGoodPMT->SetFillColor(kBlue);
  hAllGoodPMT->Draw();
  
  caaa->cd(3);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllRawBar->SetFillColor(kRed);
  hAllRawBar->Draw();
  
  caaa->cd(4);
  //gPad->SetLogy();
  //hs4->SetMinimum(0.);
  hAllGoodBar->SetFillColor(kBlue);
  hAllGoodBar->Draw();

  return caa;
}

void plotPMTRates(Int_t mymodule=1, Int_t mylayer=1, Int_t choice = 1){

  Int_t offsetl = (mylayer-1)*1344 + (mymodule-1)*224;
  Int_t offsetr = (mylayer-1)*1344 + 672 + (mymodule-1)*224;
  Int_t xcpos = 50 + (choice-1)*400;

  TCanvas *caPMT = new TCanvas(TString::Format("PMT Rates Left %d",choice), TString::Format("PMT Rates Left %d",choice), xcpos,50,400,550);
  caPMT->Divide(2,7,0.01,0.01,0);
  TCanvas *caPMTT = new TCanvas(TString::Format("PMT Rates Right %d",choice), TString::Format("PMT Rates Right %d",choice), xcpos,650,400,550);
  caPMTT->Divide(2,7,0.01,0.01,0);

  Double_t histymax = 12500;
  Double_t rate_levelu = 1800; // Expect 600 for 50000 events @ 5uA
  Double_t rate_levell = 600; // Expect 600 for 50000 events @ 5uA
  TLine *lineu = new TLine(0, rate_levelu, 2688, rate_levelu);
  TLine *linel = new TLine(0, rate_levell, 2688, rate_levell);

  for (Int_t ii=0; ii<14; ii++) {
    caPMT->cd(ii+1);
    gPad->SetLogy();
    gPad->DrawFrame(offsetl+ii*16,1.,offsetl+ii*16+16,histymax);
    hAllRawPMT->Draw("sames");

    //draw unused pixels on each PMT plot
    double xmin = offsetl + ii * 16;
    double xmax = xmin + 16;
    for (double x : missingPixelBins) {
      if (x < xmin || x >= xmax) continue;

      int bin = hAllRawPMT->FindBin(x);
      double xlow = hAllRawPMT->GetBinLowEdge(bin);
      double xup = xlow + hAllRawPMT->GetBinWidth(bin);
      double yup = hAllRawPMT->GetBinContent(bin);
      
      TBox *box = new TBox(xlow, 0, xup, yup);
      box->SetFillColor(kBlack);
      box->SetFillStyle(1001);
      box->Draw("sames");
    }
	
    lineu->Draw("sames");
    linel->Draw("sames");
  }
  caPMT->Update();

  for (Int_t ii=0; ii<14; ii++) {
    caPMTT->cd(ii+1);
    gPad->SetLogy();
    gPad->DrawFrame(offsetr+ii*16,1.,offsetr+ii*16+16.,histymax);
    hAllRawPMT->Draw("sames");

    double xmin = offsetr + ii * 16;
    double xmax = xmin + 16;
    for (double x : missingPixelBins) {
      if (x < xmin || x >= xmax) continue;

      int bin = hAllRawPMT->FindBin(x);
      double xlow = hAllRawPMT->GetBinLowEdge(bin);
      double xup = xlow + hAllRawPMT->GetBinWidth(bin);
      double yup = hAllRawPMT->GetBinContent(bin);
      
      TBox *boxx = new TBox(xlow, 0, xup, yup);
      boxx->SetFillColor(kBlack);
      boxx->SetFillStyle(1001);
      boxx->Draw("sames");
    }
    
    lineu->Draw("sames");
    linel->Draw("sames");
  }
  caPMTT->Update();
  
 
  return;


}
TCanvas *plotBarRates(){

  TCanvas *cabar = new TCanvas("BarRates", "BarRates", 50,50,800,800);
  cabar->Divide(3,4,0.01,0.01,0);

  Double_t histymax = 200000.0;
  Double_t rate_level = 600; // Expect 600 for 50000 events @ 5uA
  TLine *line = new TLine(0, rate_level, 168, rate_level);

  cabar->cd(1);
  gPad->SetLogy();
  gPad->DrawFrame(0.,1.,14.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(2);
  gPad->SetLogy();
  gPad->DrawFrame(14.,1.,28.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(3);
  gPad->SetLogy();
  gPad->DrawFrame(28.,1.,42.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(4);
  gPad->SetLogy();
  gPad->DrawFrame(42.,1.,56.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(5);
  gPad->SetLogy();
  gPad->DrawFrame(56.,1.,70.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(6);
  gPad->SetLogy();
  gPad->DrawFrame(70.,1.,84.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(7);
  gPad->SetLogy();
  gPad->DrawFrame(84.,1.,98.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(8);
  gPad->SetLogy();
  gPad->DrawFrame(98.,1.,112.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(9);
  gPad->SetLogy();
  gPad->DrawFrame(112.,1.,126.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  cabar->cd(10);
  gPad->SetLogy();
  gPad->DrawFrame(126.,1.,140.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
 
  cabar->cd(11);
  gPad->SetLogy();
  gPad->DrawFrame(140.,1.,154.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");
  
  cabar->cd(12);
  gPad->SetLogy();
  gPad->DrawFrame(154.,1.,168.,histymax);
  hAllRawBar->Draw("sames");
  line->Draw("sames");

  return cabar;


}

TCanvas *plotGoodTDC2D(){

  TCanvas *cac = new TCanvas("all2d", "all2d", 50,50,800,800);
  cac->Divide(2,2,0.01,0.01,0);
  
  cac->cd(1);
  gPad->SetLogz();
  h2AllGoodLe->Draw("colz");
  cac->cd(2);
  gPad->SetLogz();
  h2AllGoodTe->Draw("colz");
  cac->cd(3);
  gPad->SetLogz();
  h2AllGoodTot->Draw("colz");

  return cac;
}

TCanvas *plotRefTDC() {

  TCanvas *cbb = new TCanvas("ref", "ref", 850,50, 1200,800);
  cbb->Divide(2,2,0.01,0.01,0);

  cbb->cd(1);
  gPad->SetLogy();
  hRefRawLe->Draw();
  cbb->cd(2);
  gPad->SetLogy();
  hRefRawTe->Draw();
  cbb->cd(3);
  gPad->SetLogy();
  hRefRawTot->Draw();
  cbb->cd(4);
  gPad->SetLogy();
  hRefRawPMT->Draw();

  return cbb;
}


TCanvas *plotCDetTDC(){

  TCanvas *canvas[NumHalfModules];
  for (int cmodule=1;cmodule<=NumModules;cmodule++) {
   for (int layer=1;layer<=NumLayers;layer++) {
    for (int side=1;side<=NumSides;side++){
     int xposition1 = 100 + (layer-1)*1000 + (side-1)*500;
     int yposition1 = 100 + (NumModules-cmodule)*400;
     int xposition2 = 150 + (layer-1)*1000 + (side-1)*500;
     int yposition2 = 150 + (NumModules-cmodule)*400;

     int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
     int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;

     TString cname;
     cname.Form("c1_%d_%d_%d",cmodule,layer,side);
     TCanvas *c1 = new TCanvas(cname, cname, xposition1,yposition1,380,280);
     c1->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c1;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c1->cd(ii+1);
        if (ii == 1) {
		cout << "side_group_plot = " << side_group_plot << endl;
	}
	
	hRawTe[elemID_start + ii ]->Draw();
	
     }

     cname.Form("c2_%d_%d_%d",cmodule,layer,side);
     TCanvas *c2 = new TCanvas(cname, cname, xposition2,yposition2,380,280);
     c2->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
     canvas[side_group_plot] = c2;

     for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c2->cd(ii+1);
        if (ii == 1) {
                cout << "side_group_plot = " << side_group_plot << endl;
        }

        hRawLe[elemID_start + ii + 1]->Draw();

     }


    }
   }
  }

  return canvas[0];

}

TCanvas *plotHalfModule(int cmodule = 1, int side = 1, int layer = 1){


  TCanvas *c4 = new TCanvas("c4", "c4", 50,50,1000,1200);
  c4->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);
  TCanvas *c4b = new TCanvas("c4b", "c4b", 1050,50,1000,1200);
  c4b->Divide(NumPaddles,NumBars, 0.01, 0.01, 0);

  int side_group_plot = 6*(layer-1) + 2*(cmodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (cmodule-1)*NumBars*NumPaddles;


  for (int ii = 0; ii < NumPaddles*NumBars; ii++) {

        c4->cd(ii+1);
  	hRawTe[elemID_start + ii ]->Draw();
	c4b->cd(ii+1);
  	hRawLe[elemID_start + ii ]->Draw();
  }

  return c4;

}

TCanvas *plotBarTDC(int bar = 39, int side = 1, int layer = 1){

  int mymodule = (bar-1)/NumBars+1; //mymodule in this case represnts top, middle, or bottom as opposed todetector labels whe
  int paddle_start = (bar - (mymodule-1)*NumBars - 1)*NumPaddles;
  std::cout << "mymodule = " << mymodule << "  paddle_start = " << paddle_start << std::endl;

  int side_group = 6*(layer-1) + 2*(mymodule-1) + (side-1);
  int elemID_start = (layer-1)*NumModules*NumBars*NumPaddles*NumLayers + (side-1)*NumModules*NumBars*NumPaddles + (mymodule-1)*NumBars*NumPaddles + paddle_start;


  std::cout << "side_group = " << side_group << std::endl;

  TCanvas *c3 = new TCanvas("c3", "c3", 150,150,600,450);
  c3->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333 = new TCanvas("c333", "c333", 750,150,600,450);
  c333->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c3a = new TCanvas("c3a", "c3a", 150,650,600,450);
  c3a->Divide(4,4, 0.01, 0.01, 0);
  TCanvas *c333a = new TCanvas("c333a", "c333a", 750,650,600,450);
  c333a->Divide(4,4, 0.01, 0.01, 0);

  for (int ii = 0; ii < NumPaddles; ii++) {

        c3->cd(ii+1);
  	hRawLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c333->cd(ii+1);
  	hGoodLe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c3a->cd(ii+1);
  	hRawTe[elemID_start + ii ]->Draw();
  }
  for (int ii = 0; ii < NumPaddles; ii++) {

        c333a->cd(ii+1);
  	hGoodTe[elemID_start + ii ]->Draw();
  }

  return c3;

}

TCanvas *plotTOTvsLE(){

  TCanvas *c123 = new TCanvas("c123", "c123", 50,50,1000,1000);

  c123->cd();
  h2TDCTOTvsLE->Draw("colz");

  return c123;

}

TCanvas *plotTOTvsXDiff(){

  TCanvas *c1234 = new TCanvas("c1234", "c1234", 50,50,1000,1000);
  c1234->Divide(2,2,0.01,0.01,0);

  c1234->cd(1);
  h2TOTvsXDiff1->Draw("colz");
  c1234->cd(2);
  h2TOTvsXDiff2->Draw("colz");
  c1234->cd(3);
  h2LEvsXDiff1->Draw("colz");
  c1234->cd(4);
  h2LEvsXDiff2->Draw("colz");

  return c1234;

}


TCanvas *plotRowColLayer(){


  TCanvas *c5 = new TCanvas("c5", "c5", 50,50,800,800);
  c5->Divide(4,2, 0.01, 0.01, 0);

  c5->cd(1);
  gPad->SetLogy();
  hRowLayer1Side1->Draw();
  c5->cd(2);
  gPad->SetLogy();
  hRowLayer1Side2->Draw();
  c5->cd(3);
  gPad->SetLogy();
  hRowLayer2Side1->Draw();
  c5->cd(4);
  gPad->SetLogy();
  hRowLayer2Side2->Draw();
  c5->cd(5);
  gPad->SetLogy();
  hRow->Draw();
  c5->cd(6);
  gPad->SetLogy();
  hCol->Draw();
  c5->cd(7);
  gPad->SetLogy();
  hLayer->Draw();
  c5->cd(8);
  gPad->SetLogy();
  hHitPMT->Draw();

  return c5;

}

TCanvas *plotNhits(){
  TCanvas *c55 = new TCanvas("c55", "c5", 50,50,800,800);
  c55->Divide(3,4, 0.01, 0.01, 0);

  c55->cd(1);
  hnhits1->Draw();
  c55->cd(2);
  hngoodhits1->Draw();
  c55->cd(3);
  hngoodTDChits1->Draw();
  c55->cd(4);
  hnhits2->Draw();
  c55->cd(5);
  hngoodhits2->Draw();
  c55->cd(6);
  hngoodTDChits2->Draw();
  
  c55->cd(7);
  hnpaddles->Draw();
  c55->cd(8);
  hngoodpaddles->Draw();
  c55->cd(9);
  hngoodTDCpaddles->Draw();
  
  c55->cd(10);
  hnhits_ev->Draw();
  c55->cd(11);
  hngoodhits_ev->Draw();
  c55->cd(12);
  hngoodTDChits_ev->Draw();

  return c55;

}

TCanvas *plotTDC2d(){


  TCanvas *c6 = new TCanvas("c6", "c6", 50,50,800,800);
  c6->Divide(2,2, 0.01, 0.01, 0);

  c6->cd(1);
  h2d_RawLE->Draw();
  c6->cd(2);
  h2d_RawTE->Draw();
  c6->cd(3);
  h2d_RawTot->Draw();
  c6->cd(4);
  h2d_Mult->Draw();

  return c6;

}

auto *plotEECalCDet() {

  TCanvas *c1717 = new TCanvas("c1717", "c1717", 50,50,800,800);
  c1717->Divide(1,2, 0.01, 0.01, 0);

  c1717->cd(1);
  hEECalCDet1->Draw();
  c1717->cd(2);
  hEECalCDet2->Draw();

  return c1717;
}

auto plotXYZ(){

   TCanvas *c7 = new TCanvas("c7", "c7", 800,800);
   c7->Draw();
   TPad *p1 = new TPad("p1","p1",0.05,0.0,0.45,1.0);
   p1->Draw();
   p1->Divide(1,3);

   p1->cd(1);
   gPad->SetLogy();
   hHitX->Draw();

   p1->cd(2);
   gPad->SetLogy();
   hHitY->Draw();
   
   p1->cd(3);
   gPad->SetLogy();
   hHitZ->Draw();

   c7->cd(0);
   TPad *p2 = new TPad("p1","p1",0.55,0.0,0.95,1.0);
   p2->Draw();
   p2->Divide(2,1);

   p2->cd(1);
   gPad->SetLogz();
   hHitXY1->Draw("colz");
   p2->cd(2);
   gPad->SetLogz();
   hHitXY2->Draw("colz");


  return c7;

}

auto plotXYECalCDet(){

   TCanvas *c8 = new TCanvas("c8", "c7", 1200,1200);
   c8->Divide(3,3);

   c8->cd(1);
   gPad->SetLogz();
   hXECalCDet1->Draw("colz");

   c8->cd(2);
   gPad->SetLogz();
   hXECalCDet2->Draw("colz");
   
   c8->cd(3);
   hXECal->Draw();
   
   c8->cd(4);
   hYECal->Draw();
   
   c8->cd(5);
    // Define the Gaussian + constant background function
    TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

    // Set initial parameters:
    // [0] amplitude, [1] mean, [2] sigma, [3] constant background
    fitFunc->SetParameters(hXDiffECalCDet1->GetMaximum(), 0.02, 0.01, hXDiffECalCDet1->GetMinimum());
    fitFunc->SetParNames("Amplitude", "Mean", "Sigma", "Background");

    // Optional: set limits on the parameters if needed
    fitFunc->SetParLimits(1, -0.05, 0.05);  // constrain the mean near 0.02
    fitFunc->SetParLimits(2, 0.001, 0.05);  // positive sigma

    // Fit the histogram
    hXDiffECalCDet1->Fit(fitFunc, "R");  // "R" = use function range only

    // Draw the result
    hXDiffECalCDet1->Draw();
    fitFunc->Draw("same");
    // Extract Gaussian parameters
    double A = fitFunc->GetParameter(0); // Amplitude
    double mu = fitFunc->GetParameter(1); // Mean
    double sigma = fitFunc->GetParameter(2); // Sigma
    double bg = fitFunc->GetParameter(3); // Background level

    // Define integration limits
    double x_min = mu - 3*sigma;
    double x_max = mu + 3*sigma;

    // Signal: integral of the Gaussian part only over ±3σ
    TF1* gausOnly = new TF1("gausOnly", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
    gausOnly->SetParameters(A, mu, sigma);
    double signal = gausOnly->Integral(x_min, x_max);

    // Noise: integral of background over same range
    double noise = bg * (x_max - x_min);

    // Compute signal-to-noise ratio
    double snr = (noise > 0) ? signal / noise : 0;

    std::cout << "Signal (Gaussian, ±3σ): " << signal << std::endl;
    std::cout << "Noise (Background, ±3σ): " << noise << std::endl;
    std::cout << "Signal-to-Noise Ratio: " << snr << std::endl;
   

   c8->cd(6);
    // Define the Gaussian + constant background function
    TF1* fitFunc2 = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", -0.12, 0.15);

    // Set initial parameters:
    // [0] amplitude, [1] mean, [2] sigma, [3] constant background
    fitFunc2->SetParameters(hXDiffECalCDet2->GetMaximum(), 0.02, 0.01, hXDiffECalCDet2->GetMinimum());
    fitFunc2->SetParNames("Amplitude", "Mean", "Sigma", "Background");

    // Optional: set limits on the parameters if needed
    fitFunc2->SetParLimits(1, -0.05, 0.05);  // constrain the mean near 0.02
    fitFunc2->SetParLimits(2, 0.001, 0.05);  // positive sigma

    // Fit the histogram
    hXDiffECalCDet2->Fit(fitFunc2, "R");  // "R" = use function range only

    // Draw the result
    hXDiffECalCDet2->Draw();
    fitFunc2->Draw("same");

    // Extract Gaussian parameters
    double A2 = fitFunc2->GetParameter(0); // Amplitude
    double mu2 = fitFunc2->GetParameter(1); // Mean
    double sigma2 = fitFunc2->GetParameter(2); // Sigma
    double bg2 = fitFunc2->GetParameter(3); // Background level

    // Define integration limits
    double x_min2 = mu2 - 3*sigma2;
    double x_max2 = mu2 + 3*sigma2;

    // Signal: integral of the Gaussian part only over ±3σ
    TF1* gausOnly2 = new TF1("gausOnly2", "[0]*exp(-0.5*((x-[1])/[2])^2)", x_min, x_max);
    gausOnly2->SetParameters(A2, mu2, sigma2);
    double signal2 = gausOnly2->Integral(x_min2, x_max2);

    // Noise: integral of background over same range
    double noise2 = bg2 * (x_max2 - x_min2);

    // Compute signal-to-noise ratio
    double snr2 = (noise2 > 0) ? signal2 / noise2 : 0;

    std::cout << "Signal (Gaussian, ±3σ): " << signal2 << std::endl;
    std::cout << "Noise (Background, ±3σ): " << noise2 << std::endl;
    std::cout << "Signal-to-Noise Ratio: " << snr2 << std::endl;

   c8->cd(7);
   hYECalCDet1->Draw();
   
   c8->cd(8);
   hYECalCDet2->Draw();

   
   c8->cd(9);
   hXYECal->Draw("colz");
   

  return c8;
}
