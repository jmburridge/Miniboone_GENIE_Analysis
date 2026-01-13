// ============================================================================
// THIS MACRO IS FOR PROCESSING MINIBOONE DATA READY FOR THE 
// FORWARD FOLDING MACRO TO TEST THE ANALYSIS CHAIN. 

// This macro tests the response matrices and selections by passing MiniBooNE data 
// through its own response matrices. If done correctly the output should match the
// original MiniBooNE data distributions: https://arxiv.org/pdf/1805.12028
// 
// Total Truth Energy (NuMomT) Spectra created for:
//    • NuE backgrounds
//    • Pi0 backgrounds
//    • NCDelta backgrounds
//
// Uses the SAME MiniBooNE MC files and SAME background classification function:
//     StackedBkgdType_t StackHistoBkgd(...)
//     unsigned sp::Pi0Details(...)
//
// Binned in 11-bin MiniBooNE LEE energy binning.
// Writes single ROOT output file containing:
//     h_true_nue_LEE
//     h_true_pi0_LEE
//     h_true_ncdelta_LEE
//
// This output file can then be passed through the forward folding test macro and 
// then through the building macro to recereate the MiniBooNE data distributions.
//
// Output ROOT file saved to:
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Forward_Folding_macros/
// ============================================================================
#include "config.h"

#include <iostream>
#include <vector>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TStyle.h"

#include COMBINED_TYPES_H
#include COMBINED_FUNCTIONS_H
#include COMBINED_FUNCTIONS_CXX

using namespace sp;

// ----------------------------------------------------------------------------
// Pointer safety
// ----------------------------------------------------------------------------
bool CheckPointers(std::vector<int>* FSPType,
                   std::vector<float>* Vx, std::vector<float>* Vy, std::vector<float>* Vz,
                   std::vector<float>* MomX, std::vector<float>* MomY,
                   std::vector<float>* MomZ, std::vector<float>* MomT)
{
    return (FSPType && Vx && Vy && Vz && MomX && MomY && MomZ && MomT);
}

// ----------------------------------------------------------------------------
// MAIN MACRO
// ----------------------------------------------------------------------------
void process_trees()
{
    // ------------------------------------------------------------
    // LEE binning (same as in π0, Nue, NCDelta response macros)
    // ------------------------------------------------------------
    const int nbinsLEE = 11;
    double lee_bins[nbinsLEE+1] =
        {0.2, 0.3, 0.375, 0.475, 0.55, 0.675, 0.8,
         0.95, 1.1, 1.25, 1.5, 3.0};

    // ------------------------------------------------------------
    // Output histograms
    // ------------------------------------------------------------
    TH1D* h_true_nue_LEE =
        new TH1D("h_true_nue_LEE",
                 "Total True E_{#nu} (NuE backgrounds);E_{#nu} [GeV];Events",
                 nbinsLEE, lee_bins);

    TH1D* h_true_pi0_LEE =
        new TH1D("h_true_pi0_LEE",
                 "Total True E_{#nu} (Pi0 backgrounds);E_{#nu} [GeV];Events",
                 nbinsLEE, lee_bins);

    TH1D* h_true_ncdelta_LEE =
        new TH1D("h_true_ncdelta_LEE",
                 "Total True E_{#nu} (NC#Delta backgrounds);E_{#nu} [GeV];Events",
                 nbinsLEE, lee_bins);

    long long totalEntries   = 0;
    long long count_nue      = 0;
    long long count_pi0      = 0;
    long long count_ncdelta  = 0;

    // ------------------------------------------------------------
    // Loop over MiniBooNE files
    // ------------------------------------------------------------
    for (int fileIndex = 1; fileIndex <= 10; ++fileIndex)
    {
        std::stringstream ss;
        ss << OSC_MC_PREFIX << fileIndex << ".root";

        TFile* f = TFile::Open(ss.str().c_str());
        if (!f || f->IsZombie())
        {
            std::cerr << "Warning: could not open file " << ss.str() << std::endl;
            continue;
        }

        TTree* t = static_cast<TTree*>(f->Get("MiniBooNE_CCQE"));
        if (!t)
        {
            std::cerr << "Warning: TTree MiniBooNE_CCQE not found in "
                      << ss.str() << std::endl;
            f->Close();
            continue;
        }

        // --------------------------------------------------------
        // Branches
        // --------------------------------------------------------
        int   NFSP        = 0;
        int   NUANCEChan  = 0;
        int   NuType      = 0;
        int   NuParentID  = 0;
        float NuMomT      = 0.0;
        float Weight      = 0.0;
        bool  PassOsc     = false; // no selection applied, but branch is read

        std::vector<int>*   FSPType = nullptr;
        std::vector<float> *Vx = nullptr, *Vy = nullptr, *Vz = nullptr;
        std::vector<float> *MomX = nullptr, *MomY = nullptr, *MomZ = nullptr, *MomT = nullptr;

        t->SetBranchAddress("NFSP",       &NFSP);
        t->SetBranchAddress("FSPType",    &FSPType);
        t->SetBranchAddress("VertexX",    &Vx);
        t->SetBranchAddress("VertexY",    &Vy);
        t->SetBranchAddress("VertexZ",    &Vz);
        t->SetBranchAddress("MomX",       &MomX);
        t->SetBranchAddress("MomY",       &MomY);
        t->SetBranchAddress("MomZ",       &MomZ);
        t->SetBranchAddress("MomT",       &MomT);
        t->SetBranchAddress("NUANCEChan", &NUANCEChan);
        t->SetBranchAddress("NuType",     &NuType);
        t->SetBranchAddress("NuMomT",     &NuMomT);
        t->SetBranchAddress("NuParentID", &NuParentID);
        t->SetBranchAddress("Weight",     &Weight);
        t->SetBranchAddress("PassOsc",    &PassOsc);

        const int N = t->GetEntries();
        totalEntries += N;

        // --------------------------------------------------------
        // Event loop
        // --------------------------------------------------------
        for (int i = 0; i < N; ++i)
        {
            t->GetEntry(i);

            if (!CheckPointers(FSPType, Vx, Vy, Vz, MomX, MomY, MomZ, MomT))
                continue;

            // ----------------------------------------------------
            // Determine if event has a true pi0 using Pi0Details
            // ----------------------------------------------------
            unsigned npi0 = sp::Pi0Details(NFSP,
                                           *FSPType,
                                           *Vx, *Vy, *Vz,
                                           *MomX, *MomY, *MomZ, *MomT);

            bool Event_is_pi0 = (npi0 > 0);

            // ----------------------------------------------------
            // Classify background type using MiniBooNE logic
            // ----------------------------------------------------
            StackedBkgdType_t bkg =
                StackHistoBkgd(false,                    // Event_is_dirt
                               Event_is_pi0,             // Event_is_pi0 from truth
                               static_cast<NuanceType_t>(NUANCEChan),
                               static_cast<NuType_t>(NuType),
                               static_cast<GEANT3Type_t>(NuParentID));

            // ----------------------------------------------------
            // Fill appropriate histogram by classification
            // ----------------------------------------------------
            if (bkg == kBKGD_PI0)
            {
                h_true_pi0_LEE->Fill(NuMomT, Weight);
                ++count_pi0;
            }
            else if (bkg == kBKGD_DELTA)
            {
                h_true_ncdelta_LEE->Fill(NuMomT, Weight);
                ++count_ncdelta;
            }
            else if (bkg == kBKGD_NUEPIP ||
                     bkg == kBKGD_NUEKP  ||
                     bkg == kBKGD_NUEK0)
            {
                h_true_nue_LEE->Fill(NuMomT, Weight);
                ++count_nue;
            }
        }

        f->Close();
    }

    // ------------------------------------------------------------
    // Output ROOT file
    // ------------------------------------------------------------
    std::string outfile = std::string(TEST_CACHE_DIR) + "/" + PATH_FF_TREES_NUANCE;
    TFile fout(outfile.c_str(), "RECREATE");

    h_true_nue_LEE->Write();
    h_true_pi0_LEE->Write();
    h_true_ncdelta_LEE->Write();

    fout.Close();

    // ------------------------------------------------------------
    // Summary
    // ------------------------------------------------------------
    std::cout << "=============================================\n";
    std::cout << "  LEE Total Truth Background Spectra Written\n";
    std::cout << "  Total events processed : " << totalEntries   << "\n";
    std::cout << "  NuE     events counted : " << count_nue      << "\n";
    std::cout << "  Pi0     events counted : " << count_pi0      << "\n";
    std::cout << "  NCDelta events counted : " << count_ncdelta  << "\n";
    std::cout << "  Output file: " << outfile                   << "\n";
    std::cout << "=============================================\n";
}

