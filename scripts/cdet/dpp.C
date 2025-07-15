#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>

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
        return;
    }
    h->SetTitle("heep.dpp with cuts;#delta p/p;Counts");
    h->Draw();

    //TFile *outFile = new TFile("heep_dpp_all_files.root", "RECREATE");
    //h->Write();
    //outFile->Close();

    printf("Done!");
}
