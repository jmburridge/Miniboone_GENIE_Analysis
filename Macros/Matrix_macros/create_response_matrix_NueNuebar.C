// ============================================================================
// Final Clean νe Response Matrix Macro (Improved)
// Energy-only, MicroBooNE-style full response matrix
//
// νe signal defined via StackHistoBkgd categories:
//   kBKGD_NUEPIP, kBKGD_NUEKP, kBKGD_NUEK0
//
// Standardised ROOT outputs (matching π0 macro style):
//   - response_nue.root       (13-bin "technote" νe binning)
//   - response_nue_LEE.root   (11-bin MiniBooNE LEE binning)
//
// ROOT files saved to:
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Nue/Nue_Root_Files/
//
// PNGs saved to:
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Nue/Nue_Histograms/
//
// Styling parallels π0 macro:
//   - kInvertedDarkBodyRadiator palette
//   - Stat boxes on 1D
//   - Improved margins
// ============================================================================

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"

#include "../../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedTypes.h"
#include "../../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.h"
#include "../../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.cxx"

using namespace sp;

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------
bool CheckPointers(std::vector<int>* FSPType, std::vector<float>* Vx,
                   std::vector<float>* Vy, std::vector<float>* Vz,
                   std::vector<float>* MomX, std::vector<float>* MomY,
                   std::vector<float>* MomZ, std::vector<float>* MomT)
{
    return (FSPType && Vx && Vy && Vz && MomX && MomY && MomZ && MomT);
}

void CopyAndLabelIndexed(TH2D* src, TH2D* idx)
{
    int nx = src->GetNbinsX();
    int ny = src->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
            idx->SetBinContent(ix, iy, src->GetBinContent(ix, iy));

    for (int i = 1; i <= nx; ++i)
        idx->GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
    for (int j = 1; j <= ny; ++j)
        idx->GetYaxis()->SetBinLabel(j, std::to_string(j).c_str());
}

// -----------------------------------------------------------------------------
// Main macro
// -----------------------------------------------------------------------------
void create_response_matrix_NueNuebar()
{
    // -------------------------------------------------------------------------
    // Global style
    // -------------------------------------------------------------------------
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    gStyle->SetOptStat(1110);      // entries, mean, RMS
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatX(0.88);
    gStyle->SetStatY(0.88);

    // -------------------------------------------------------------------------
    // Binning: 13-bin νe scheme + 11-bin LEE scheme
    // -------------------------------------------------------------------------
    const int nbins = 13;

    double reco_bins[nbins + 1] = {
        0.200, 0.300, 0.375, 0.475, 0.550, 0.675, 0.800,
        0.950, 1.100, 1.250, 1.500, 1.750, 2.000, 2.500
    };

    double true_bins[nbins + 1] = {
        0.200, 0.250, 0.300, 0.350, 0.400,
        0.450, 0.500, 0.600, 0.800, 1.000,
        1.500, 2.000, 2.500, 3.000
    };

    const int lee_nbins = 11;
    double lee_bins[lee_nbins + 1] = {
        0.2, 0.3, 0.375, 0.475, 0.55,
        0.675, 0.8, 0.95, 1.1, 1.3,
        1.5, 3.0
    };

    // -------------------------------------------------------------------------
    // Standard-binning histograms (13-bin)
    // -------------------------------------------------------------------------
    TH1D* h_total_true = new TH1D("h_total_true",
        "Total True E_{#nu_{e}};E_{#nu_{e}} [GeV];Events",
        nbins, true_bins);

    TH1D* h_pass_true = new TH1D("h_pass_true",
        "Passed True E_{#nu_{e}};E_{#nu_{e}} [GeV];Events",
        nbins, true_bins);

    TH1D* h_pass_reco = new TH1D("h_pass_reco",
        "Passed Reco E_{#nu_{e}};Reco E_{#nu_{e}} [GeV];Events",
        nbins, reco_bins);

    TH2D* h_response_E = new TH2D("h_response_E",
        "Response Matrix (Reco vs True E_{#nu_{e}});Reco E_{#nu_{e}} [GeV];True E_{#nu_{e}} [GeV]",
        nbins, reco_bins,
        nbins, true_bins);

    TH1D* h_eff_nue = new TH1D("h_eff_nue",
        "ν_{e} Efficiency;True E_{#nu_{e}} [GeV];Efficiency",
        nbins, true_bins);

    // -------------------------------------------------------------------------
    // LEE-binned histograms (11-bin)
    // -------------------------------------------------------------------------
    TH1D* h_total_true_LEE = new TH1D("h_total_true_LEE",
        "Total True E_{#nu_{e}} (LEE);E_{#nu_{e}} [GeV];Events",
        lee_nbins, lee_bins);

    TH1D* h_pass_true_LEE = new TH1D("h_pass_true_LEE",
        "Passed True E_{#nu_{e}} (LEE);E_{#nu_{e}} [GeV];Events",
        lee_nbins, lee_bins);

    TH2D* h_response_E_LEE = new TH2D("h_response_E_LEE",
        "LEE Response Matrix (Reco vs True E_{#nu_{e}});Reco E_{#nu_{e}} [GeV];True E_{#nu_{e}} [GeV]",
        lee_nbins, lee_bins,
        lee_nbins, lee_bins);

    // -------------------------------------------------------------------------
    // Event counters
    // -------------------------------------------------------------------------
    long long totalEntries = 0;
    long long nue_count    = 0;

    // -------------------------------------------------------------------------
    // Loop over input files
    // -------------------------------------------------------------------------
    for (int fileIndex = 1; fileIndex <= 10; ++fileIndex)
    {
        std::stringstream ss;
        ss << "../../MiniBooNEDatasets2023/output_osc_mc_detail_"
           << fileIndex << ".root";

        TFile* f = TFile::Open(ss.str().c_str());
        if (!f || f->IsZombie()) continue;

        TTree* t = (TTree*)f->Get("MiniBooNE_CCQE");
        if (!t) { f->Close(); continue; }

        int   NFSP, NUANCEChan, NuType, NuParentID;
        float Energy, RecoEnuQE, Weight, NuMomT;
        bool  PassOsc;

        std::vector<int>*   FSPType = nullptr;
        std::vector<float>* Vx      = nullptr;
        std::vector<float>* Vy      = nullptr;
        std::vector<float>* Vz      = nullptr;
        std::vector<float>* MomX    = nullptr;
        std::vector<float>* MomY    = nullptr;
        std::vector<float>* MomZ    = nullptr;
        std::vector<float>* MomT    = nullptr;

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
        t->SetBranchAddress("Energy",     &Energy);
        t->SetBranchAddress("RecoEnuQE",  &RecoEnuQE);
        t->SetBranchAddress("Weight",     &Weight);
        t->SetBranchAddress("PassOsc",    &PassOsc);

        int N = t->GetEntries();
        totalEntries += N;

        for (int i = 0; i < N; ++i)
        {
            t->GetEntry(i);
            if (!CheckPointers(FSPType, Vx, Vy, Vz, MomX, MomY, MomZ, MomT))
                continue;

            // νe classification using StackHistoBkgd
            StackedBkgdType_t bkgd_type =
                StackHistoBkgd(false,                  // Event_is_dirt
                               false,                  // Event_is_pi0
                               (NuanceType_t)NUANCEChan,
                               (NuType_t)NuType,
                               (GEANT3Type_t)NuParentID);

            bool isNueSignal =
                (bkgd_type == kBKGD_NUEPIP ||
                 bkgd_type == kBKGD_NUEKP  ||
                 bkgd_type == kBKGD_NUEK0);

            if (!isNueSignal)
                continue;

            ++nue_count;

            // Fill standard-binning truth spectrum
            h_total_true->Fill(NuMomT, Weight);
            // Fill LEE-binning truth spectrum
            h_total_true_LEE->Fill(NuMomT, Weight);

            if (PassOsc) {
                // Standard-binning "passed" truth & reco
                h_pass_true->Fill(NuMomT, Weight);
                h_pass_reco->Fill(RecoEnuQE, Weight);
                h_response_E->Fill(RecoEnuQE, NuMomT, Weight);

                // LEE-binning
                h_pass_true_LEE->Fill(NuMomT, Weight);
                h_response_E_LEE->Fill(RecoEnuQE, NuMomT, Weight);
            }
        }

        f->Close();
    }

    // -------------------------------------------------------------------------
    // Normalise standard 13×13 response matrix (column by column in true bin)
    // -------------------------------------------------------------------------
    for (int iy = 1; iy <= nbins; ++iy) {
        double den = h_total_true->GetBinContent(iy);
        if (den <= 0.0) continue;

        for (int ix = 1; ix <= nbins; ++ix) {
            double val = h_response_E->GetBinContent(ix, iy);
            h_response_E->SetBinContent(ix, iy, val / den);
        }
    }

    // -------------------------------------------------------------------------
    // Normalise LEE 11×11 response matrix
    // -------------------------------------------------------------------------
    for (int iy = 1; iy <= lee_nbins; ++iy) {
        double den = h_total_true_LEE->GetBinContent(iy);
        if (den <= 0.0) continue;

        for (int ix = 1; ix <= lee_nbins; ++ix) {
            double val = h_response_E_LEE->GetBinContent(ix, iy);
            h_response_E_LEE->SetBinContent(ix, iy, val / den);
        }
    }

    // -------------------------------------------------------------------------
    // νe efficiency vs true energy (standard 13-bin)
    // -------------------------------------------------------------------------
    for (int i = 1; i <= nbins; ++i) {
        double num = h_pass_true->GetBinContent(i);
        double den = h_total_true->GetBinContent(i);

        if (den > 0.0) {
            double eff = num / den;
            double err = std::sqrt(eff * (1.0 - eff) / den);
            h_eff_nue->SetBinContent(i, eff);
            h_eff_nue->SetBinError(i, err);
        }
    }

    // -------------------------------------------------------------------------
    // Indexed response matrices
    // -------------------------------------------------------------------------
    TH2D* h_response_E_idx =
        new TH2D("h_response_E_idx", ";Reco E bin;True E bin",
                 nbins, 0, nbins, nbins, 0, nbins);
    CopyAndLabelIndexed(h_response_E, h_response_E_idx);

    TH2D* h_response_E_idx_LEE =
        new TH2D("h_response_E_idx_LEE", ";Reco E bin;True E bin",
                 lee_nbins, 0, lee_nbins, lee_nbins, 0, lee_nbins);
    CopyAndLabelIndexed(h_response_E_LEE, h_response_E_idx_LEE);

    // -------------------------------------------------------------------------
    // Output paths
    // -------------------------------------------------------------------------
    const char* png_dir =
        "/exp/uboone/app/users/jburridg/Geometry/Analysis/Nue/Nue_Histograms/";
    const char* root_dir =
        "/exp/uboone/app/users/jburridg/Geometry/Analysis/Nue/Nue_Root_Files/";

    // -------------------------------------------------------------------------
    // Plotting utilities
    // -------------------------------------------------------------------------
    auto SaveMatrix = [&](TH2* h, const char* fname)
    {
        TCanvas* c = new TCanvas(fname, fname, 900, 800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.08);
        h->LabelsOption("h");
        h->Draw("COLZ");
        std::string out = std::string(png_dir) + fname + ".png";
        c->SaveAs(out.c_str());
    };

    auto Save1D = [&](TH1* h, const char* fname)
    {
        TCanvas* c = new TCanvas(fname, fname, 800, 800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.15);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.08);
        h->Draw("HIST E");
        std::string out = std::string(png_dir) + fname + ".png";
        c->SaveAs(out.c_str());
    };

    // -------------------------------------------------------------------------
    // Save PNGs
    // -------------------------------------------------------------------------
    // Standard
    SaveMatrix(h_response_E,       "Nue_ResponseMatrix_Energy");
    SaveMatrix(h_response_E_idx,   "Nue_ResponseMatrix_Energy_Indexed");

    // LEE
    SaveMatrix(h_response_E_LEE,     "Nue_ResponseMatrix_Energy_LEE");
    SaveMatrix(h_response_E_idx_LEE, "Nue_ResponseMatrix_Energy_LEE_Indexed");

    // 1D
    Save1D(h_total_true,      "Nue_TrueEnergy_Total");
    Save1D(h_pass_true,       "Nue_TrueEnergy_Passed");
    Save1D(h_pass_reco,       "Nue_RecoEnergy_Passed");
    Save1D(h_eff_nue,         "Nue_Efficiency_TrueEnergy");
    Save1D(h_total_true_LEE,  "Nue_TrueEnergy_Total_LEE");
    Save1D(h_pass_true_LEE,   "Nue_TrueEnergy_Passed_LEE");

    // -------------------------------------------------------------------------
    // Standardised export: 13-bin νe (response_nue.root)
    // -------------------------------------------------------------------------
    {
        TH1D* h_true_export = (TH1D*)h_total_true->Clone("h_true");
        TH1D* h_pass_export = (TH1D*)h_pass_true->Clone("h_pass");
        TH1D* h_reco_export = (TH1D*)h_pass_reco->Clone("h_reco");

        // Following π0 macro convention: smear_unorm is cloned from normalised matrix
        TH2D* h_smear_unorm = (TH2D*)h_response_E->Clone("h_smear_unorm");
        TH2D* h_resp        = (TH2D*)h_response_E->Clone("h_resp");

        TH2D* h_resp_index =
            new TH2D("h_resp_index", "Response indexed",
                     nbins, 0, nbins, nbins, 0, nbins);
        CopyAndLabelIndexed(h_resp, h_resp_index);

        std::string outfile = std::string(root_dir) + "response_nue.root";
        TFile out(outfile.c_str(), "RECREATE");
        h_true_export->Write();
        h_pass_export->Write();
        h_reco_export->Write();
        h_smear_unorm->Write();
        h_resp->Write();
        h_resp_index->Write();
        out.Close();
    }

    // -------------------------------------------------------------------------
    // Standardised export: LEE-binned νe (response_nue_LEE.root)
    // -------------------------------------------------------------------------
    {
        TH1D* h_true_LEE_export = (TH1D*)h_total_true_LEE->Clone("h_true_LEE");
        TH1D* h_pass_LEE_export = (TH1D*)h_pass_true_LEE->Clone("h_pass_LEE");

        // Build reco LEE distribution from response (sum over true bins)
        TH1D* h_reco_LEE =
            new TH1D("h_reco_LEE", "Reco Energy (LEE-binned)", lee_nbins, lee_bins);

        for (int ix = 1; ix <= lee_nbins; ++ix) {
            double sum = 0.0;
            for (int iy = 1; iy <= lee_nbins; ++iy)
                sum += h_response_E_LEE->GetBinContent(ix, iy);
            h_reco_LEE->SetBinContent(ix, sum);
        }

        TH2D* h_smear_LEE_unorm = (TH2D*)h_response_E_LEE->Clone("h_smear_LEE_unorm");
        TH2D* h_resp_LEE        = (TH2D*)h_response_E_LEE->Clone("h_resp_LEE");

        TH2D* h_resp_index_LEE =
            new TH2D("h_resp_index_LEE", "LEE Response indexed",
                     lee_nbins, 0, lee_nbins, lee_nbins, 0, lee_nbins);
        CopyAndLabelIndexed(h_resp_LEE, h_resp_index_LEE);

        std::string outfile = std::string(root_dir) + "response_nue_LEE.root";
        TFile out(outfile.c_str(), "RECREATE");
        h_true_LEE_export->Write();
        h_pass_LEE_export->Write();
        h_reco_LEE->Write();
        h_smear_LEE_unorm->Write();
        h_resp_LEE->Write();
        h_resp_index_LEE->Write();
        out.Close();
    }

    // -------------------------------------------------------------------------
    // Summary
    // -------------------------------------------------------------------------
    std::cout
        << "[NUE] Final νe response matrices (standard + LEE) completed.\n"
        << "      Entries processed: " << totalEntries << "\n"
        << "      νe events kept:    " << nue_count   << std::endl;
}
