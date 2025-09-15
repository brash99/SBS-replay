#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>

// Scalars (1D vectors)
std::vector<double> vheep_dpp;
std::vector<double> vheep_dt_ADC;
std::vector<double> vheep_ecalo;
std::vector<double> vheep_eprime_eth;
std::vector<double> vheep_dxECAL;
std::vector<double> vearm_ecal_x;
std::vector<int>    vsbs_gemFPP_ntrack;
std::vector<double> vheep_dyECAL;

// Arrays (2D vectors)
std::vector<std::vector<double>> vsbs_tr_vz;
std::vector<std::vector<double>> vsbs_gemFPP_sclose;
std::vector<std::vector<int>>    vsbs_gemFT_nhits;
std::vector<std::vector<int>>    vsbs_gemFT_ngoodhits;

template<typename T>
std::vector<T> fill2D(const TTreeReaderArray<T>& arr) {
  std::vector<T> tmp;
  int n = arr.GetSize();
  tmp.reserve(n);
  for (int i=0; i<n; i++) tmp.push_back(arr[i]);
  return tmp;
}


/* Main Routine */
void dpp(const char* dirname = ".") {
    TChain *chain = new TChain("T");

    // Open directory
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    if (!files) {
        printf("No files found in directory %s\n", dirname);
        return;
    }

    // Loop over all files
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file = (TSystemFile*)next())) {
        fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
            TString fullpath = TString(dirname) + "/" + fname;
            printf("Adding file: %s\n", fullpath.Data());
            chain->Add(fullpath);
        }
    }

    const char *cut = "abs(heep.dt_ADC)<10 && abs(sbs.tr.vz[0]+0.1)<0.18 && heep.ecalo/heep.eprime_eth > 0.7 && abs(heep.dxECAL - 0.01 + 0.025 * earm.ecal.x) < 0.05 && sbs.gemFPP.track.ntrack > 0 && abs(heep.dyECAL - 0.01) < 0.06 && sbs.gemFPP.track.sclose[0] < 0.01 && (sbs.gemFT.track.nhits[0] > 4 || sbs.gemFT.track.ngoodhits[0] > 2)";

    chain->Draw("heep.dpp >> hHeepDpp(100, -0.1, 0.1)", cut, "hist");

    TH1F *h = (TH1F*)gDirectory->Get("hHeepDpp");
    if (!h) {
        printf("Histogram hHeepDpp not found!\n");
    }
    /*
    // Set up TTreeReader
    TTreeReader reader(chain);

    // Scalars
    TTreeReaderValue<double> heep_dpp(reader, "heep.dpp");
    TTreeReaderValue<double> heep_dt_ADC(reader, "heep.dt_ADC");
    TTreeReaderValue<double> heep_ecalo(reader, "heep.ecalo");
    TTreeReaderValue<double> heep_eprime_eth(reader, "heep.eprime_eth");
    TTreeReaderValue<double> heep_dxECAL(reader, "heep.dxECAL");
    TTreeReaderValue<double> earm_ecal_x(reader, "earm.ecal.x");
    TTreeReaderValue<int> sbs_gemFPP_ntrack(reader, "sbs.gemFPP.track.ntrack");
    TTreeReaderValue<double> heep_dyECAL(reader, "heep.dyECAL");

    // Arrays
    TTreeReaderArray<double> sbs_tr_vz(reader, "sbs.tr.vz");
    TTreeReaderArray<double> sbs_gemFPP_sclose(reader, "sbs.gemFPP.track.sclose");
    TTreeReaderArray<int> sbs_gemFT_nhits(reader, "sbs.gemFT.track.nhits");
    TTreeReaderArray<int> sbs_gemFT_ngoodhits(reader, "sbs.gemFT.track.ngoodhits");

    // Fill vector from TChain
    while (reader.Next()) {
        // Scalars — push_back the dereferenced values
        vheep_dpp.push_back(*heep_dpp);
        vheep_dt_ADC.push_back(*heep_dt_ADC);
        vheep_ecalo.push_back(*heep_ecalo);
        vheep_eprime_eth.push_back(*heep_eprime_eth);
        vheep_dxECAL.push_back(*heep_dxECAL);
        vearm_ecal_x.push_back(*earm_ecal_x);
        vsbs_gemFPP_ntrack.push_back(*sbs_gemFPP_ntrack);
        vheep_dyECAL.push_back(*heep_dyECAL);

        //2D arrays
        vsbs_tr_vz.push_back(fill2D(sbs_tr_vz));
        vsbs_gemFPP_sclose.push_back(fill2D(sbs_gemFPP_sclose));
        vsbs_gemFT_nhits.push_back(fill2D(sbs_gemFT_nhits));
        vsbs_gemFT_ngoodhits.push_back(fill2D(sbs_gemFT_ngoodhits));
    }
        */

    std::cout << "Filled vheep_dpp with " << vheep_dpp.size() << " entries.\n";

    const char *cut = "abs(heep.dt_ADC)<10 && abs(sbs.tr.vz[0]+0.1)<0.18 && heep.ecalo/heep.eprime_eth > 0.7 && abs(heep.dxECAL - 0.01 + 0.025 * earm.ecal.x) < 0.05 && sbs.gemFPP.track.ntrack > 0 && abs(heep.dyECAL - 0.01) < 0.06 && sbs.gemFPP.track.sclose[0] < 0.01 && (sbs.gemFT.track.nhits[0] > 4 || sbs.gemFT.track.ngoodhits[0] > 2)";

    chain->Draw("heep.dpp >> hheep_dpp(100, -0.1, 0.1)", cut, "hist");

    TH1F *h = (TH1F*)gDirectory->Get("hheep_dpp");
    if (!h) {
        printf("Histogram hheep_dpp not found!\n");
        return;
    }
    h->SetTitle("heep.dpp with cuts;#delta p/p;Counts");
    h->Draw();

    //TFile *outFile = new TFile("heep_dpp_all_files.root", "RECREATE");
    //h->Write();
    //outFile->Close();

    printf("Done!");
}


