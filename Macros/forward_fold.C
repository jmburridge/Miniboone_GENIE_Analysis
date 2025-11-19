///////////////////////// INFORMATION ///////////////////////////
//
//A macro to forward fold GENIE MC truth spectra through MiniBooNE response matrices
//
//This macro accepts 2 input .root files: 
//		1. A .root file containing all response matrices for the
//		   background categories.
//		2. A .root file containing all the MC truth spectra for 
//		   NC Delta, Pi0, and Nue background categories 
//
//This macro creates 1 output: 
//		1. A .root file containing the Reconstructed spectra for 
//		   NC Delta, Pi0, and Nue background categories. 
//
//
///////////////////////////////////////////////////////////////

static const char* PATH_BACKGROUND_TREES =
    "../DigitisePlot/background_event_trees.root";

static const char* PATH_PI0_MATRIX =
    "../Matrix_macros/response_matrix_root_files/response_pi0_LEE.root";

static const char* PATH_DELTA_MATRIX =
    "../Matrix_macros/response_matrix_root_files/response_ncdeltas_LEE.root";

static const char* PATH_NUE_MATRIX =
   "../Matrix_macros/response_matrix_root_files/response_nues_LEE.root";

// ------------------------------------------------------------
// Response histogram names — EDIT IF NEEDED
// ------------------------------------------------------------
static const char* RESP_PI0_NAME   = "h_resp_LEE";
static const char* RESP_DELTA_NAME = "h_resp_LEE";
static const char* RESP_NUE_NAME   = "h_resp_LEE";

// ------------------------------------------------------------
// Output file (also hardcoded)
// ------------------------------------------------------------
static const char* PATH_OUTPUT =
    "forward_folds_output.root";

// ------------------------------------------------------------
// LEE binning scheme (truth-energy variable binning)
// ------------------------------------------------------------
static const int NBINS = 11;
double LEE_bins[NBINS+1] = {
    0.2, 0.3, 0.375, 0.475, 0.55,
    0.675, 0.8, 0.95, 1.1, 1.3,
    1.5, 3.0
};

// ------------------------------------------------------------
// PassOsc requirement
// ------------------------------------------------------------
static const bool REQUIRE_PASSOSC = true;

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <cmath>

// ------------------------------------------------------------
// Build truth histogram with LEE bins
// ------------------------------------------------------------
TH1D* BuildTruthLEE(TTree* t, const char* hname)
{
    if (!t) {
        std::cerr << "[BuildTruthLEE] ERROR: null tree " << hname << std::endl;
        return nullptr;
    }

    Float_t NuMomT = 0;
    Float_t Weight = 1;
    Bool_t  PassOsc = true;

    t->SetBranchStatus("*",0);
    t->SetBranchStatus("NuMomT",   1);
    t->SetBranchStatus("Weight",  1);
    t->SetBranchStatus("PassOsc", 1);

    t->SetBranchAddress("NuMomT",&NuMomT);
    t->SetBranchAddress("Weight",&Weight);
    t->SetBranchAddress("PassOsc",&PassOsc);

    TH1D* h = new TH1D(hname, hname, NBINS, LEE_bins);
    h->Sumw2();
    h->Reset();

    Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; i++) {
        t->GetEntry(i);
        if (REQUIRE_PASSOSC && !PassOsc) continue;
        h->Fill(NuMomT, Weight);
    }

    return h;
}

// ------------------------------------------------------------
// Forward-folding: reco_i = sum_j R(i,j) * truth_j
// Assumes X=reco, Y=true
// ------------------------------------------------------------
TH1D* ForwardFold(const TH1D* hTruth, const TH2D* hResp, const char* outName)
{
    if (!hTruth || !hResp) {
        std::cerr << "[ForwardFold] ERROR: null inputs." << std::endl;
        return nullptr;
    }

    int nTruth = hResp->GetNbinsY();
    int nReco  = hResp->GetNbinsX();

    double* recoEdges = new double[nReco+1];
    for (int b=1; b<=nReco+1; ++b)
        recoEdges[b-1] = hResp->GetXaxis()->GetBinLowEdge(b);

    TH1D* hReco = new TH1D(outName, outName, nReco, recoEdges);
    hReco->Sumw2();
    hReco->Reset();
    delete [] recoEdges;

    for (int j=1; j<=nTruth; j++) {
        double Tj = hTruth->GetBinContent(j);
        if (Tj == 0) continue;

        for (int i=1; i<=nReco; i++) {
            double R = hResp->GetBinContent(i,j);
            double contrib = Tj * R;

            double old = hReco->GetBinContent(i);
            double err = hReco->GetBinError(i);

            hReco->SetBinContent(i, old + contrib);
            hReco->SetBinError(i, std::sqrt(err*err + std::fabs(contrib)));
        }
    }

    return hReco;
}

// ------------------------------------------------------------
// Main function (HARD-CODED PATHS)
// ------------------------------------------------------------
void ForwardFoldMiniBooNE()
{
    // ------------------------
    // Load background trees
    // ------------------------
    TFile* fT = TFile::Open(PATH_BACKGROUND_TREES, "READ");
    if (!fT || fT->IsZombie()) {
        std::cerr << "ERROR: Cannot open background tree file: "
                  << PATH_BACKGROUND_TREES << std::endl;
        return;
    }

    TTree* tPi0    = (TTree*)fT->Get("pi0");
    TTree* tDelta  = (TTree*)fT->Get("delta");
    TTree* tNuePip = (TTree*)fT->Get("nuepip");
    TTree* tNueKp  = (TTree*)fT->Get("nuekp");
    TTree* tNueK0  = (TTree*)fT->Get("nuek0");

    // ------------------------
    // Build truth histograms
    // ------------------------
    TH1D* hPi0   = BuildTruthLEE(tPi0, "hTruth_pi0");
    TH1D* hDelta = BuildTruthLEE(tDelta, "hTruth_delta");

    TH1D* hNuepip = BuildTruthLEE(tNuePip, "hTruth_nuepip");
    TH1D* hNuekp  = BuildTruthLEE(tNueKp, "hTruth_nuekp");
    TH1D* hNuek0  = BuildTruthLEE(tNueK0, "hTruth_nuek0");

    // Sum intrinsic nu_e components
    TH1D* hNue = nullptr;
    if (hNuepip) {
        hNue = (TH1D*)hNuepip->Clone("hTruth_nue");
        if (hNuekp) hNue->Add(hNuekp);
        if (hNuek0) hNue->Add(hNuek0);
    }

    // ------------------------
    // Load response matrices
    // ------------------------
    TFile* fRpi0   = TFile::Open(PATH_PI0_MATRIX, "READ");
    TFile* fRdelta = TFile::Open(PATH_DELTA_MATRIX, "READ");
    TFile* fRnue   = TFile::Open(PATH_NUE_MATRIX, "READ");

    TH2D* Rpi0   = (TH2D*)fRpi0  ->Get(RESP_PI0_NAME);
    TH2D* Rdelta = (TH2D*)fRdelta->Get(RESP_DELTA_NAME);
    TH2D* Rnue   = (TH2D*)fRnue  ->Get(RESP_NUE_NAME);

    // ------------------------
    // Forward folding
    // ------------------------
    TH1D* FFpi0   = (hPi0   && Rpi0)   ? ForwardFold(hPi0,   Rpi0,   "hReco_pi0")   : nullptr;
    TH1D* FFdelta = (hDelta && Rdelta) ? ForwardFold(hDelta, Rdelta, "hReco_delta") : nullptr;
    TH1D* FFnue   = (hNue   && Rnue)   ? ForwardFold(hNue,   Rnue,   "hReco_nue")   : nullptr;

    // ------------------------
    // Write output
    // ------------------------
    TFile* fOut = TFile::Open(PATH_OUTPUT, "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "ERROR: Cannot create output file: " << PATH_OUTPUT << std::endl;
        return;
    }

    if (hPi0)   hPi0->Write();
    if (hDelta) hDelta->Write();
    if (hNuepip) hNuepip->Write();
    if (hNuekp)  hNuekp->Write();
    if (hNuek0)  hNuek0->Write();
    if (hNue)    hNue->Write();

    if (FFpi0)   FFpi0->Write();
    if (FFdelta) FFdelta->Write();
    if (FFnue)   FFnue->Write();

    fOut->Close();

    std::cout << "[ForwardFold_LEE] Finished. Output written to:\n"
              << PATH_OUTPUT << std::endl;
}
