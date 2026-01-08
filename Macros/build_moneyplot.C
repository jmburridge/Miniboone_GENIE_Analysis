////////////////////// INFORMATION /////////////////////
//
//A macro to stack the reconstructed  spectra into a
//'MiniBooNE LEE' - style plot.
//
//
///////////////////////////////////////////////////////
// ============================================================================
// MiniBooNE Background Stack Constructor (Corrected Final Version)
// - Dirt & Other: weighted + POT scaled
// - Pi0 / Delta / intrinsic nue: from forward-folding, POT scaled (check this)
// - All MC histograms rebinned into identical 200–3000 MeV bins
// - MiniBooNE 'dirt' and 'other' hists from text files (MiniBooNEDatasets2023), converted to Events/MeV
// - Stacks all baclground categories correctly
// ============================================================================

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"

#include "../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.h"
#include "../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.cxx"
#include "../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedTypes.h"

#include <fstream>
#include <iostream>
#include <vector>

// ============================================================================
// BINNING
// ============================================================================
static const int NBINS = 11;
double LEE_bins[NBINS+1] = {
    200, 300, 375, 475, 550,
    675, 800, 950, 1100, 1250,
    1500, 3000
};

// ============================================================================
// Divide by bin width
// ============================================================================
void DivideByBinWidth(TH1D* h)
{
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        double c = h->GetBinContent(i);
        double e = h->GetBinError(i);
        double w = h->GetBinWidth(i);
        if (w > 0) {
            h->SetBinContent(i, c / w);
            h->SetBinError(i,   e / w);
        }
    }
}

// ============================================================================
// Read MiniBooNE reference (Events per bin)
// ============================================================================
TH1D* Read_MB_Text(const char* path)
{
    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::cerr << "ERROR opening MB text file\n";
        return nullptr;
    }

    TH1D* h = new TH1D("MiniBooNE_Total", "MiniBooNE_Total", NBINS, LEE_bins);
    h->Sumw2();

    double v; int b = 1;
    while (fin >> v && b <= NBINS) {
        h->SetBinContent(b, v);
        h->SetBinError(b, std::sqrt(std::max(v, 0.0)));
        b++;
    }

    DivideByBinWidth(h);
    return h;
}

// ============================================================================
// Build hist from TTree (Dirt/Other)
// ============================================================================
TH1D* Build_From_TTree(TTree* t, const char* name)
{
    if (!t) return nullptr;
    Float_t RecoE, Weight;

    t->SetBranchAddress("RecoE", &RecoE);
    t->SetBranchAddress("Weight", &Weight);

    TH1D* h = new TH1D(name, name, NBINS, LEE_bins);
    h->Sumw2();

    Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; i++) {
        t->GetEntry(i);
        h->Fill(RecoE * 1000.0, Weight); // GeV → MeV
    }
    return h;
}

// ============================================================================
// Load forward-folded histogram
// ============================================================================
TH1D* LoadFF(TFile* f, const char* name)
{
    TH1D* h = (TH1D*)f->Get(name);
    if (!h) {
        std::cerr << "ERROR: missing FF histogram " << name << "\n";
        return nullptr;
    }
    TH1D* c = (TH1D*)h->Clone((TString)name + "_clone");
    c->SetDirectory(nullptr);
    c->Sumw2();
    return c;
}

// ============================================================================
// *** Correct Fix ***
//  Rebin ANY forward-folded histogram into the LEE bin array
//  converting x-axis from GeV → MeV
// ============================================================================
TH1D* RebinToLEE(const TH1D* src, const char* newName)
{
    TH1D* dst = new TH1D(newName, newName, NBINS, LEE_bins);
    dst->Sumw2();

    for (int i = 1; i <= src->GetNbinsX(); i++) {
        double xGeV = src->GetBinCenter(i);
        double xMeV = xGeV * 1000.0;

        double val = src->GetBinContent(i);
        double err = src->GetBinError(i);

        int b = dst->FindBin(xMeV);
        if (b < 1 || b > NBINS) continue;

        dst->AddBinContent(b, val);
        double eold = dst->GetBinError(b);
        dst->SetBinError(b, std::sqrt(eold*eold + err*err));
    }
    return dst;
}

// ============================================================================
// External Dirt
// ============================================================================
void Add_External_Dirt(TH1D* hDirt)
{
    TFile* df = TFile::Open("../MiniBooNEDatasets2023/miniboone_mc_dirt_postccqe.root");
    if (!df || df->IsZombie()) return;

    TTree* t = (TTree*)df->Get("MiniBooNBE_CCQE");
    if (!t) t = (TTree*)df->Get("MiniBooNE_CCQE");
    if (!t) return;

    int NFSP, NUANCEChan, NuType, NuParentID;
    float RecoEnuQE, Weight;
    bool PassOsc;

    std::vector<int>* FSPType_v = nullptr;
    std::vector<float>* X = nullptr; std::vector<float>* Y = nullptr;
    std::vector<float>* Z = nullptr; std::vector<float>* MX = nullptr;
    std::vector<float>* MY = nullptr; std::vector<float>* MZ = nullptr;
    std::vector<float>* MT = nullptr;

    t->SetBranchAddress("NFSP", &NFSP);
    t->SetBranchAddress("FSPType", &FSPType_v);
    t->SetBranchAddress("VertexX", &X);
    t->SetBranchAddress("VertexY", &Y);
    t->SetBranchAddress("VertexZ", &Z);
    t->SetBranchAddress("MomX", &MX);
    t->SetBranchAddress("MomY", &MY);
    t->SetBranchAddress("MomZ", &MZ);
    t->SetBranchAddress("MomT", &MT);
    t->SetBranchAddress("NUANCEChan", &NUANCEChan);
    t->SetBranchAddress("NuType", &NuType);
    t->SetBranchAddress("NuParentID", &NuParentID);
    t->SetBranchAddress("RecoEnuQE", &RecoEnuQE);
    t->SetBranchAddress("Weight", &Weight);
    t->SetBranchAddress("PassOsc", &PassOsc);

    Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; i++) {
        t->GetEntry(i);
        if (!PassOsc) continue;

        unsigned isPi0 =
            sp::Pi0Details(NFSP, *FSPType_v, *X, *Y, *Z, *MX, *MY, *MZ, *MT);

        auto bkg =
            sp::StackHistoBkgd(true, isPi0,
                               (sp::NuanceType_t)NUANCEChan,
                               (sp::NuType_t)NuType,
                               (sp::GEANT3Type_t)NuParentID);

        if (bkg == sp::kBKGD_DIRT)
            hDirt->Fill(RecoEnuQE * 1000.0, Weight);
    }
    df->Close();
}

// ============================================================================
// MAIN
// ============================================================================
void build_miniboone_plot()
{
    gStyle->SetOptStat(0);

    // Load Dirt/Other TTrees
    TFile* fTrees = TFile::Open("processed_trees_for_FF.root");
    TH1D* hDirt  = Build_From_TTree((TTree*)fTrees->Get("dirt"),  "hDirt");
    TH1D* hOther = Build_From_TTree((TTree*)fTrees->Get("other"), "hOther");

    //Add_External_Dirt(hDirt);

    // Load forward-folds (raw)
    TFile* fFF = TFile::Open("forward_folds_output.root");
    TH1D* hPi0_raw   = LoadFF(fFF, "hReco_pi0");
    TH1D* hDelta_raw = LoadFF(fFF, "hReco_delta");
    TH1D* hNue_raw   = LoadFF(fFF, "hReco_nue");

    // Rebin FF into correct binning
    TH1D* hPi0   = RebinToLEE(hPi0_raw,   "hPi0");
    TH1D* hDelta = RebinToLEE(hDelta_raw, "hDelta");
    TH1D* hNue   = RebinToLEE(hNue_raw,   "hNue");

    // Load MiniBooNE reference
    TH1D* hMB = Read_MB_Text("../Archive/DigitisePlot/nu_lee2018_numode_background.txt");

    // Convert MC → Events/MeV
    for (TH1D* h : {hDirt, hOther, hDelta, hNue, hPi0}) //hDelta, hNue, hPi0
        DivideByBinWidth(h);

    // Scale EVERYTHING by the same factor as Dirt/Other
    double scale = 12.84 / 41.1;
    hDirt->Scale(scale);
    hOther->Scale(scale);

    // New: apply same scaling to FF categories
    hPi0->Scale(scale);
    hDelta->Scale(scale);
    hNue->Scale(scale);


    // Colors
    hDirt->SetFillColor(kGray+1);     hDirt->SetLineColor(kBlack);
    hOther->SetFillColor(kAzure-9);   hOther->SetLineColor(kBlue+2);
    hPi0->SetFillColor(kOrange-3);    hPi0->SetLineColor(kOrange+7);
    hDelta->SetFillColor(kRed-9);     hDelta->SetLineColor(kRed+1);
    hNue->SetFillColor(kGreen-7);     hNue->SetLineColor(kGreen+3);

    // Stack
    THStack* hs = new THStack(
        "stack",
        "MiniBooNE Background Stack;Reconstructed Neutrino Energy (MeV);Events/MeV");
    
    hs->SetMaximum(6);
    hs->Add(hOther);
    hs->Add(hDirt);
    hs->Add(hDelta);
    hs->Add(hPi0);
    hs->Add(hNue);
   
   

    // Canvas
    TCanvas* c = new TCanvas("c","LEE stack",900,700);

    hs->Draw("HIST");
    hs->SetMinimum(0);
    hs->SetMaximum(5);
    hs->Draw("HIST");

    // MiniBooNE curve
    hMB->SetLineColor(kBlack);
    hMB->SetLineWidth(3);
    hMB->SetLineStyle(2);   // dashed line
    hMB->SetMarkerStyle(0); // no markers
    hMB->SetMaximum(6);
    hMB->Draw("HIST SAME"); // draw as line, not points


    // Legend
    TLegend* leg = new TLegend(0.60,0.55,0.88,0.88);
    leg->AddEntry(hDirt,  "Dirt", "f");
    leg->AddEntry(hOther, "Other", "f");
    leg->AddEntry(hPi0,   "#pi^{0}", "f");
    leg->AddEntry(hDelta, "NC #Delta", "f");
    leg->AddEntry(hNue,   "Intrinsic #nu_{e}", "f");
    leg->AddEntry(hMB, "Expected background", "l");
    leg->Draw();

    // Save
    c->SaveAs("Stacked_LEE_backgrounds_with_MB.png");

    // ROOT output
    TFile fout("Stacked_LEE_backgrounds_with_MB.root","RECREATE");
    hDirt->Write(); hOther->Write(); hPi0->Write();
    hDelta->Write(); hNue->Write(); hMB->Write();
    hs->Write();
    fout.Close();
}
