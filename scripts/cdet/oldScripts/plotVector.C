// Refactored version of PlotHVScanHighCurrent.C
// Key changes:
// 1. Vectors collect data in the event loop
// 2. Histograms are declared and filled only right before plotting
// 3. Plotting routines (e.g. plotAllTDC) handle filling and drawing

#include <TROOT.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>

// Constants
const int nTdc = 2704;

// Event data vectors
std::vector<double> goodLeVec, goodTeVec, goodTotVec;
std::vector<int> goodIDVec;
std::vector<double> goodXVec, goodYVec, goodZVec;
std::vector<double> ecalXVec, ecalYVec, ecalEVec;
std::vector<int> layerVec, sideVec, rowVec, colVec;

// Placeholder histogram pointers
TH1F *hGoodLe = nullptr;
TH1F *hGoodTe = nullptr;
TH1F *hGoodTot = nullptr;
std::vector<TH1F*> hGoodTotBarVec;
TH2F *hXYECal = nullptr;
TH1F *hXECal = nullptr, *hYECal = nullptr, *hEECal = nullptr;
TH1F *hHitX = nullptr, *hHitY = nullptr, *hHitZ = nullptr;
TH1F *hLayer = nullptr, *hCol = nullptr, *hRow = nullptr;

void FillHistograms_GoodTDC() {
    hGoodLe = new TH1F("hGoodLe", "Good LE Times", 100, 0, 60);
    hGoodTe = new TH1F("hGoodTe", "Good TE Times", 100, 0, 60);
    hGoodTot = new TH1F("hGoodTot", "Good TOT", 100, 2, 45);
    for (size_t i = 0; i < goodLeVec.size(); ++i) {
        hGoodLe->Fill(goodLeVec[i]);
        hGoodTe->Fill(goodTeVec[i]);
        hGoodTot->Fill(goodTotVec[i]);
    }
}

void FillHistograms_TotByBar() {
    hGoodTotBarVec.resize(nTdc);
    for (int i = 0; i < nTdc; ++i) {
        hGoodTotBarVec[i] = new TH1F(Form("hGoodTot_Bar%d", i), Form("TOT Bar %d", i), 100, 2, 45);
    }
    for (size_t i = 0; i < goodTotVec.size(); ++i) {
        int id = goodIDVec[i];
        if (id >= 0 && id < nTdc) {
            hGoodTotBarVec[id]->Fill(goodTotVec[i]);
        }
    }
}

void FillHistograms_Positions() {
    hHitX = new TH1F("hHitX", "X Hit", 200, -2, 2);
    hHitY = new TH1F("hHitY", "Y Hit", 200, -0.5, 0.5);
    hHitZ = new TH1F("hHitZ", "Z Hit", 200, 7.5, 8.5);
    for (size_t i = 0; i < goodXVec.size(); ++i) {
        hHitX->Fill(goodXVec[i]);
        hHitY->Fill(goodYVec[i]);
        hHitZ->Fill(goodZVec[i]);
    }
}

void FillHistograms_ECal() {
    hXYECal = new TH2F("hXYECal", "ECal XY", 200, -2.0, 2.0, 200, -2.0, 2.0);
    hXECal = new TH1F("hXECal", "ECal X", 200, -1.5, 1.5);
    hYECal = new TH1F("hYECal", "ECal Y", 200, -1.0, 1.0);
    hEECal = new TH1F("hEECal", "ECal Energy", 200, 0.0, 20.0);
    for (size_t i = 0; i < ecalXVec.size(); ++i) {
        hXYECal->Fill(ecalYVec[i], ecalXVec[i]);
        hXECal->Fill(ecalXVec[i]);
        hYECal->Fill(ecalYVec[i]);
        hEECal->Fill(ecalEVec[i]);
    }
}

void FillHistograms_LayerColRow() {
    hLayer = new TH1F("hLayer", "Layer Number", 3, 0, 3);
    hCol = new TH1F("hCol", "Side Number", 3, 0, 3);
    hRow = new TH1F("hRow", "Row Number", 700, 0, 700);
    for (size_t i = 0; i < layerVec.size(); ++i) {
        hLayer->Fill(layerVec[i]);
        hCol->Fill(sideVec[i]);
        hRow->Fill(rowVec[i]);
    }
}

void PlotHVScanHighCurrent(Int_t RunNumber1=3562, Int_t nevents=50000) {
    TChain *T = new TChain("T");
    TString infile = Form("/path/to/cdet_%d.root", RunNumber1);
    T->Add(infile);

    double le, te, tot, x, y, z, ex, ey, ee;
    int id, lyr, col, row, side;

    T->SetBranchAddress("tdc_le", &le);
    T->SetBranchAddress("tdc_te", &te);
    T->SetBranchAddress("tdc_tot", &tot);
    T->SetBranchAddress("tdc_id", &id);
    T->SetBranchAddress("x", &x);
    T->SetBranchAddress("y", &y);
    T->SetBranchAddress("z", &z);
    T->SetBranchAddress("ecal_x", &ex);
    T->SetBranchAddress("ecal_y", &ey);
    T->SetBranchAddress("ecal_e", &ee);
    T->SetBranchAddress("layer", &lyr);
    T->SetBranchAddress("col", &col);
    T->SetBranchAddress("row", &row);
    T->SetBranchAddress("side", &side);

    for (int i = 0; i < nevents; ++i) {
        T->GetEntry(i);
        goodLeVec.push_back(le);
        goodTeVec.push_back(te);
        goodTotVec.push_back(tot);
        goodIDVec.push_back(id);
        goodXVec.push_back(x);
        goodYVec.push_back(y);
        goodZVec.push_back(z);
        ecalXVec.push_back(ex);
        ecalYVec.push_back(ey);
        ecalEVec.push_back(ee);
        layerVec.push_back(lyr);
        colVec.push_back(col);
        rowVec.push_back(row);
        sideVec.push_back(side);
    }

    //plotAllTDC();
    plotTOTByBar();

    //FillHistograms_Positions();
    //FillHistograms_ECal();
    //FillHistograms_LayerColRow();

    TFile *fout = new TFile("output_plots.root", "RECREATE");
    if (hGoodLe) hGoodLe->Write();
    if (hGoodTe) hGoodTe->Write();
    if (hGoodTot) hGoodTot->Write();
    if (hHitX) hHitX->Write();
    if (hHitY) hHitY->Write();
    if (hHitZ) hHitZ->Write();
    if (hXYECal) hXYECal->Write();
    if (hXECal) hXECal->Write();
    if (hYECal) hYECal->Write();
    if (hEECal) hEECal->Write();
    if (hLayer) hLayer->Write();
    if (hCol) hCol->Write();
    if (hRow) hRow->Write();
    for (auto h : hGoodTotBarVec) if (h) h->Write();
    fout->Close();
}

void plotAllTDC() {
    FillHistograms_GoodTDC();
    TCanvas *c1 = new TCanvas("c1", "Good LE", 800, 600);
    hGoodLe->Draw();

    TCanvas *c2 = new TCanvas("c2", "Good TE", 800, 600);
    hGoodTe->Draw();

    TCanvas *c3 = new TCanvas("c3", "Good TOT", 800, 600);
    hGoodTot->Draw();
}

void plotTOTByBar() {
    FillHistograms_TotByBar();
    TCanvas *c4 = new TCanvas("c4", "TOT by Bar", 800, 600);
    for (int i = 0; i < nTdc; ++i) {
        if (hGoodTotBarVec[i]->GetEntries() > 0) {
            hGoodTotBarVec[i]->Draw("same");
        }
    }
}
