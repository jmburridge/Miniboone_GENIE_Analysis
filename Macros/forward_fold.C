///////////////////////// INFORMATION ///////////////////////////////////////////////////
//                                                                                     //
//A macro to forward fold GENIE MC truth spectra through MiniBooNE response matrices   //
//                                                                                     //
//This macro accepts 2 input .root files:                                              //
//		1. A .root file containing all response matrices for the                       //
//		   background categories.                                                      //
//		2. A .root file containing all the MC truth spectra for                        //
//		   NC Delta, Pi0, and Nue background categories                                //
//                                                                                     //
//This macro creates 1 output:                                                         //
//		1. A .root file containing the Reconstructed spectra for                       //
//		   NC Delta, Pi0, and Nue background categories.                               //
//                                                                                     //
//                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////

// ============================================================================
// ForwardFold_LEE.C  (Histogram-based version compatible with process_trees)
// ============================================================================
//
// Inputs from process_trees():
//   h_true_pi0_LEE
//   h_true_ncdelta_LEE
//   h_true_nue_LEE
//
// Categories folded:
//   pi0, delta, nue
//
// Dirt and Other are not produced by process_trees and therefore not folded.
//
// Produces: forward_folds_output.root
// ============================================================================

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <cmath>

// ============================
// File locations
// ============================
static const char* PATH_TREES =
    "truth_trees_for_FF.root";

static const char* PATH_PI0_MATRIX   =
    "../Pi0/Pi0_Root_Files/response_pi0_LEE.root";

static const char* PATH_DELTA_MATRIX =
    "../NCDelta/NCDelta_Root_Files/response_ncdelta_LEE.root";

static const char* PATH_NUE_MATRIX   =
    "../Nue/Nue_Root_Files/response_nue_LEE.root";

static const char* RESP_NAME = "h_resp_LEE";

static const char* PATH_OUTPUT =
    "forward_folds_output.root";

// ============================================================================
// Forward-fold matrix product: reco_i = Σ_j R(i,j)·truth_j
// ============================================================================
TH1D* ForwardFold(const TH1D* hTruth, const TH2D* R, const char* name)
{
    if (!hTruth || !R) return nullptr;

    const int nReco  = R->GetNbinsX();
    const int nTruth = R->GetNbinsY();

    TH1D* hReco =
        new TH1D(name, name, nReco, R->GetXaxis()->GetXbins()->GetArray());
    hReco->Sumw2();

    for (int j=1; j<=nTruth; j++) {
        const double Tj = hTruth->GetBinContent(j);
        if (Tj == 0) continue;

        for (int i=1; i<=nReco; i++) {
            const double Rij = R->GetBinContent(i, j);
            const double contrib = Tj * Rij;

            const double old    = hReco->GetBinContent(i);
            const double oldErr = hReco->GetBinError(i);

            hReco->SetBinContent(i, old + contrib);
            hReco->SetBinError(i, std::sqrt(oldErr*oldErr + std::fabs(contrib)));
        }
    }

    return hReco;
}

// ============================================================================
// MAIN
// ============================================================================
void ForwardFoldMiniBooNE()
{
    // ------------------------------------
    // Load the truth histogram file
    // ------------------------------------
    TFile* fT = TFile::Open(PATH_TREES, "READ");
    if (!fT || fT->IsZombie()) {
        std::cerr << "ERROR: cannot open " << PATH_TREES << "\n";
        return;
    }

    // --------------------------------------------------------
    // Read the histograms produced by process_trees()
    // --------------------------------------------------------
    TH1D* hPi0   = (TH1D*)fT->Get("h_true_pi0_LEE");
    TH1D* hDelta = (TH1D*)fT->Get("h_true_ncdelta_LEE");
    TH1D* hNue   = (TH1D*)fT->Get("h_true_nue_LEE");

    if (!hPi0 || !hDelta || !hNue) {
        std::cerr << "ERROR: missing one or more truth histograms\n";
        return;
    }

    // ------------------------------------
    // Load response matrices
    // ------------------------------------
    TFile* fMpi0   = TFile::Open(PATH_PI0_MATRIX, "READ");
    TFile* fMdelta = TFile::Open(PATH_DELTA_MATRIX, "READ");
    TFile* fMnue   = TFile::Open(PATH_NUE_MATRIX, "READ");

    TH2D* Rpi0   = (TH2D*)fMpi0  ->Get(RESP_NAME);
    TH2D* Rdelta = (TH2D*)fMdelta->Get(RESP_NAME);
    TH2D* Rnue   = (TH2D*)fMnue  ->Get(RESP_NAME);

    if (!Rpi0 || !Rdelta || !Rnue) {
        std::cerr << "ERROR: missing one or more response matrices\n";
        return;
    }

    // ------------------------------------
    // Forward-fold all categories
    // ------------------------------------
    TH1D* FFpi0   = ForwardFold(hPi0,   Rpi0,   "hReco_pi0");
    TH1D* FFdelta = ForwardFold(hDelta, Rdelta, "hReco_delta");
    TH1D* FFnue   = ForwardFold(hNue,   Rnue,   "hReco_nue");

    // ------------------------------------
    // Output file
    // ------------------------------------
    TFile* fout = new TFile(PATH_OUTPUT, "RECREATE");

    // Write truth histograms
    hPi0  ->Write();
    hDelta->Write();
    hNue  ->Write();

    // Write folded spectra
    FFpi0  ->Write();
    FFdelta->Write();
    FFnue  ->Write();

    fout->Close();

    std::cout << "[ForwardFoldMiniBooNE] DONE\n";
    std::cout << "Wrote: " << PATH_OUTPUT << std::endl;
}

