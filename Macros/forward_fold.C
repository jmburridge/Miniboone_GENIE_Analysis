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


#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <iostream>
#include <vector>
#include <cmath>

//=======================================================================
// User must supply correct MiniBooNE LEE bin edges here.
//=======================================================================

static const int NBINS = 20;   // Example; adjust to actual value

double LEE_bins[NBINS+1] = {
    0.0, 0.1, 0.2, 0.3, 0.4,   // Replace with official MiniBooNE LEE edges
    0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 1.1, 1.2, 1.3, 1.4,
    1.5, 1.6, 1.7, 1.8, 1.9,
    2.0
};

//=======================================================================
// Utility: compare two arrays of double bin edges
//=======================================================================
bool CheckBinEdges(const TAxis* ax, const double* refEdges, int nBins)
{
    if (ax->GetNbins() != nBins) return false;

    for (int i = 0; i <= nBins; i++) {
        double diff = std::fabs(ax->GetBinLowEdge(i+1) - refEdges[i]);
        if (diff > 1e-9) return false;
    }
    return true;
}

//=======================================================================
// Utility: Validate that a response matrix matches canonical LEE binning
//=======================================================================
bool ValidateResponseMatrix(TH2* h, const char* name)
{
    if (!h) {
        std::cerr << "ERROR: Missing response matrix: " << name << "\n";
        return false;
    }

    // Must be square
    if (h->GetNbinsX() != h->GetNbinsY()) {
        std::cerr << "ERROR: Response matrix " << name
                  << " is not square: "
                  << h->GetNbinsX() << "x" << h->GetNbinsY() << "\n";
        return false;
    }

    // Must match LEE reco binning (x-axis)
    if (!CheckBinEdges(h->GetXaxis(), LEE_bins, NBINS)) {
        std::cerr << "ERROR: Response matrix " << name
                  << " has incorrect reconstructed-energy binning.\n";
        return false;
    }

    // Must match LEE true binning (y-axis)
    if (!CheckBinEdges(h->GetYaxis(), LEE_bins, NBINS)) {
        std::cerr << "ERROR: Response matrix " << name
                  << " has incorrect true-energy binning.\n";
        return false;
    }

    std::cout << "Validated response matrix: " << name << "\n";
    return true;
}

//=======================================================================
// Utility: Validate that a truth spectrum matches canonical LEE binning
//=======================================================================
bool ValidateTruthSpectrum(TH1* h, const char* name)
{
    if (!h) {
        std::cerr << "ERROR: Missing truth spectrum: " << name << "\n";
        return false;
    }

    if (h->GetNbinsX() != NBINS) {
        std::cerr << "ERROR: Truth spectrum " << name
                  << " has wrong number of bins: "
                  << h->GetNbinsX() << " (expected " << NBINS << ")\n";
        return false;
    }

    if (!CheckBinEdges(h->GetXaxis(), LEE_bins, NBINS)) {
        std::cerr << "ERROR: Truth spectrum " << name
                  << " has incorrect true-energy binning.\n";
        return false;
    }

    std::cout << "Validated truth spectrum: " << name << "\n";
    return true;
}

//=======================================================================
// MAIN MACRO
//=======================================================================
void ForwardFoldMiniBooNE(const char* responseFile,
                          const char* truthFile,
                          const char* recoFile)
{
    std::cout << "----------------------------------------------------\n";
    std::cout << " Forward-Folding with Strict MiniBooNE LEE Binning\n";
    std::cout << "----------------------------------------------------\n";

    // Open input files
    TFile* fResp = TFile::Open(responseFile, "READ");
    if (!fResp || fResp->IsZombie()) {
        std::cerr << "ERROR: Could not open response file.\n";
        return;
    }

    TFile* fTruth = TFile::Open(truthFile, "READ");
    if (!fTruth || fTruth->IsZombie()) {
        std::cerr << "ERROR: Could not open truth spectra file.\n";
        return;
    }

    //--------------------------------------------------------------
    // Load response matrices
    //--------------------------------------------------------------
    TH2D* rm_NCDelta = (TH2D*)fResp->Get("ResponseMatrix_NCDelta");
    TH2D* rm_Pi0     = (TH2D*)fResp->Get("ResponseMatrix_Pi0");
    TH2D* rm_Nue     = (TH2D*)fResp->Get("ResponseMatrix_Nue");

    //--------------------------------------------------------------
    // Validate response matrices
    //--------------------------------------------------------------
    if (!ValidateResponseMatrix(rm_NCDelta, "ResponseMatrix_NCDelta")) return;
    if (!ValidateResponseMatrix(rm_Pi0,     "ResponseMatrix_Pi0"))     return;
    if (!ValidateResponseMatrix(rm_Nue,     "ResponseMatrix_Nue"))     return;

    //--------------------------------------------------------------
    // Load truth spectra
    //--------------------------------------------------------------
    TH1D* true_NCDelta = (TH1D*)fTruth->Get("Truth_NCDelta");
    TH1D* true_Pi0     = (TH1D*)fTruth->Get("Truth_Pi0");
    TH1D* true_Nue     = (TH1D*)fTruth->Get("Truth_Nue");

    //--------------------------------------------------------------
    // Validate truth spectra
    //--------------------------------------------------------------
    if (!ValidateTruthSpectrum(true_NCDelta, "Truth_NCDelta")) return;
    if (!ValidateTruthSpectrum(true_Pi0,     "Truth_Pi0"))     return;
    if (!ValidateTruthSpectrum(true_Nue,     "Truth_Nue"))     return;

    //--------------------------------------------------------------
    // Forward folding function: M × T
    //--------------------------------------------------------------
    auto ForwardFold = [&](TH2D* response, TH1D* truth, const char* name) {

        TH1D* reco = new TH1D(name, name, NBINS, LEE_bins);

        for (int i = 1; i <= NBINS; i++) {
            double sum = 0.0;
            for (int j = 1; j <= NBINS; j++) {
                sum += response->GetBinContent(i, j) *
                       truth->GetBinContent(j);
            }
            reco->SetBinContent(i, sum);
        }

        return reco;
    };

    //--------------------------------------------------------------
    // Perform forward folding
    //--------------------------------------------------------------
    TH1D* reco_NCDelta = ForwardFold(rm_NCDelta, true_NCDelta, "Reco_NCDelta");
    TH1D* reco_Pi0     = ForwardFold(rm_Pi0,     true_Pi0,     "Reco_Pi0");
    TH1D* reco_Nue     = ForwardFold(rm_Nue,     true_Nue,     "Reco_Nue");

    //--------------------------------------------------------------
    // Write outputs
    //--------------------------------------------------------------
    TFile* fOut = TFile::Open(recoFile, "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "ERROR: Could not create output file.\n";
        return;
    }

    reco_NCDelta->Write();
    reco_Pi0->Write();
    reco_Nue->Write();

    fOut->Close();
    fResp->Close();
    fTruth->Close();

    std::cout << "All reconstructed spectra successfully written to: "
              << recoFile << "\n";
    std::cout << "----------------------------------------------------\n";
}
