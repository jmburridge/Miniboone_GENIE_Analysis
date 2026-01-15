// ============================================================================
// This is a macro to create all matrices (Nue, NCDelta, Pi0) in one go.
// It uses MC files from MiniBooNE datasets 2023.
//
// Creating one macro to do all three to ensure orthogonal selections.
// The theory issue is that when matrices are made separately, events can be
// double-counted across channels (e.g. both Pi0 and NCDelta). By imposing
// orthogonal selections, the normalisations should be more accurate and hence
// give more accurate response matrices.
//
// ROOT files saved to:
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Pi0/Pi0_Root_Files/
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/NCDelta/NCDelta_Root_Files/
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Nue/Nue_Root_Files/
// PNGs saved to:
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Pi0/Pi0_Histograms/
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/NCDelta/NCDelta_Histograms/
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Nue/Nue_Histograms/
//
// THIS MACRO USES THE MiniBooNE 2023 MC DATASET FILES:
//   output_osc_mc_detail_1.root ... output_osc_mc_detail_10.root
// ============================================================================
#include "config.h"

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

#include COMBINED_FUNCTIONS_H
#include COMBINED_TYPES_H
#include COMBINED_FUNCTIONS_CXX

using namespace sp;

// -----------------------------------------------------------------------------
// Helpers
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

    for (int ix = 1; ix <= nx; ix++)
        for (int iy = 1; iy <= ny; iy++)
            idx->SetBinContent(ix, iy, src->GetBinContent(ix, iy));

    for (int i = 1; i <= nx; i++)
        idx->GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
    for (int j = 1; j <= ny; j++)
        idx->GetYaxis()->SetBinLabel(j, std::to_string(j).c_str());
}

// Divide bin contents by bin width → events per MeV
auto NormaliseByBinWidth = [](TH1D* h) {
    int nb = h->GetNbinsX();
    for (int i = 1; i <= nb; i++) {
        double c = h->GetBinContent(i);
        double w = h->GetBinWidth(i);
        if (w > 0) h->SetBinContent(i, c / w);
    }
    h->GetYaxis()->SetTitle("Events/MeV");
};

// -----------------------------------------------------------------------------
// Main macro
// -----------------------------------------------------------------------------
void create_all_matrices()
{
    // GLOBAL STYLE
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    gStyle->SetOptStat(1110);
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatX(0.88);
    gStyle->SetStatY(0.88);

    // -------------------------------------------------------------------------
    // Binning
    // -------------------------------------------------------------------------
    const int nbins_true_nue     = 13;
    const int nbins_reco_nue     = 13;
    const int nbins_true_ncdelta = 7;
    const int nbins_reco_ncdelta = 8;
    const int nbins_true_pi0     = 9;
    const int nbins_reco_pi0     = 9;

    // Nue binning
    double recoE_bins_nue[nbins_reco_nue+1] =
        {0.200, 0.300, 0.375, 0.475, 0.550, 0.675, 0.800,
         0.950, 1.100, 1.250, 1.500, 1.750, 2.000, 2.500};

    double trueE_bins_nue[nbins_true_nue+1] =
        {0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500,
         0.600, 0.800, 1.000, 1.500, 2.000, 2.500, 3.000};

    // NCDelta binning
    double recoE_bins_ncdelta[nbins_reco_ncdelta+1] =
        {0.200, 0.300, 0.385, 0.485, 0.550, 0.690, 0.800,
         1.000, 1.200};

    double trueE_bins_ncdelta[nbins_true_ncdelta+1] =
        {0.250, 0.500, 0.750, 1.000, 1.250, 1.500, 2.000, 3.000};

    // Pi0 binning
    double recoE_bins_pi0[nbins_reco_pi0+1] =
        {0.150, 0.200, 0.300, 0.375, 0.475, 0.550, 0.675,
         0.800, 1.000, 1.200};

    double trueE_bins_pi0[nbins_true_pi0+1] =
        {0.250, 0.500, 0.600, 0.700, 0.900, 1.100, 1.250,
         1.500, 2.0, 3.0};

    // LEE binning
    const int lee_nbins = 11;
    double lee_bins[lee_nbins+1] =
        {0.2, 0.3, 0.375, 0.475, 0.55, 0.675, 0.8,
         0.95, 1.1, 1.25, 1.5, 3.0};

    // -------------------------------------------------------------------------
    // Energy-only histograms
    // -------------------------------------------------------------------------
    TH1D* h_total_true_nue     = new TH1D("h_total_true_nue",
        "Total True E_{#nu};E_{#nu} [GeV];Events", nbins_true_nue,     trueE_bins_nue);
    TH1D* h_total_true_ncdelta = new TH1D("h_total_true_ncdelta",
        "Total True E_{#nu};E_{#nu} [GeV];Events", nbins_true_ncdelta, trueE_bins_ncdelta);
    TH1D* h_total_true_pi0     = new TH1D("h_total_true_pi0",
        "Total True E_{#nu};E_{#nu} [GeV];Events", nbins_true_pi0,     trueE_bins_pi0);

    TH1D* h_pass_true_nue     = new TH1D("h_pass_true_nue",
        "Passed True E_{#nu};E_{#nu} [GeV];Events", nbins_true_nue,     trueE_bins_nue);
    TH1D* h_pass_true_ncdelta = new TH1D("h_pass_true_ncdelta",
        "Passed True E_{#nu};E_{#nu} [GeV];Events", nbins_true_ncdelta, trueE_bins_ncdelta);
    TH1D* h_pass_true_pi0     = new TH1D("h_pass_true_pi0",
        "Passed True E_{#nu};E_{#nu} [GeV];Events", nbins_true_pi0,     trueE_bins_pi0);

    TH1D* h_pass_reco_nue     = new TH1D("h_pass_reco_nue",
        "Passed Reco E_{#nu}^{QE};Reco E_{#nu}^{QE} [GeV];Events",
        nbins_reco_nue,     recoE_bins_nue);
    TH1D* h_pass_reco_ncdelta = new TH1D("h_pass_reco_ncdelta",
        "Passed Reco E_{#nu}^{QE};Reco E_{#nu}^{QE} [GeV];Events",
        nbins_reco_ncdelta, recoE_bins_ncdelta);
    TH1D* h_pass_reco_pi0     = new TH1D("h_pass_reco_pi0",
        "Passed Reco E_{#nu}^{QE};Reco E_{#nu}^{QE} [GeV];Events",
        nbins_reco_pi0,     recoE_bins_pi0);

    TH2D* h_response_E_nue = new TH2D("h_response_E_nue",
        "Response Matrix (Reco vs True E);Reco E_{#nu}^{QE};True E_{#nu}",
        nbins_reco_nue, recoE_bins_nue, nbins_true_nue, trueE_bins_nue);

    TH2D* h_response_E_ncdelta = new TH2D("h_response_E_ncdelta",
        "Response Matrix (Reco vs True E);Reco E_{#nu}^{QE};True E_{#nu}",
        nbins_reco_ncdelta, recoE_bins_ncdelta, nbins_true_ncdelta, trueE_bins_ncdelta);

    TH2D* h_response_pi0 = new TH2D("h_response_E_pi0",
        "Response Matrix (Reco vs True E);Reco E_{#nu}^{QE};True E_{#nu}",
        nbins_reco_pi0, recoE_bins_pi0, nbins_true_pi0, trueE_bins_pi0);

    // (efficiency histos defined but not yet used)
    TH1D* h_eff_nue     = new TH1D("h_eff_nue",
        "Nue Efficiency;True E_{#nu} [GeV];Efficiency",
        nbins_true_nue,     trueE_bins_nue);
    TH1D* h_eff_ncdelta = new TH1D("h_eff_ncdelta",
        "NCDelta Efficiency;True E_{#nu} [GeV];Efficiency",
        nbins_true_ncdelta, trueE_bins_ncdelta);
    TH1D* h_eff_pi0     = new TH1D("h_eff_pi0",
        "Pi0 Efficiency;True E_{#nu} [GeV];Efficiency",
        nbins_true_pi0,     trueE_bins_pi0);

    // -------------------------------------------------------------------------
    // LEE-binned histograms
    // -------------------------------------------------------------------------
    TH1D* h_total_true_LEE_nue     = new TH1D("h_total_true_LEE_nue",
        "Total True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);
    TH1D* h_total_true_LEE_ncdelta = new TH1D("h_total_true_LEE_ncdelta",
        "Total True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);
    TH1D* h_total_true_LEE_pi0     = new TH1D("h_total_true_LEE_pi0",
        "Total True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);

    TH1D* h_pass_true_LEE_nue     = new TH1D("h_pass_true_LEE_nue",
        "Passed True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);
    TH1D* h_pass_true_LEE_ncdelta = new TH1D("h_pass_true_LEE_ncdelta",
        "Passed True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);
    TH1D* h_pass_true_LEE_pi0     = new TH1D("h_pass_true_LEE_pi0",
        "Passed True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);

    TH1D* h_pass_reco_LEE_nue     = new TH1D("h_pass_reco_LEE_nue",
        "Passed Reco E_{#nu}^{QE} (LEE);Reco E_{#nu}^{QE} [GeV];Events",
        lee_nbins, lee_bins);
    TH1D* h_pass_reco_LEE_ncdelta = new TH1D("h_pass_reco_LEE_ncdelta",
        "Passed Reco E_{#nu}^{QE} (LEE);Reco E_{#nu}^{QE} [GeV];Events",
        lee_nbins, lee_bins);
    TH1D* h_pass_reco_LEE_pi0     = new TH1D("h_pass_reco_LEE_pi0",
        "Passed Reco E_{#nu}^{QE} (LEE);Reco E_{#nu}^{QE} [GeV];Events",
        lee_nbins, lee_bins);

    TH2D* h_response_E_LEE_nue = new TH2D("h_response_E_LEE_nue",
        "LEE Response Matrix;Reco E_{#nu}^{QE};True E_{#nu}",
        lee_nbins, lee_bins, lee_nbins, lee_bins);

    TH2D* h_response_E_LEE_ncdelta = new TH2D("h_response_E_LEE_ncdelta",
        "LEE Response Matrix;Reco E_{#nu}^{QE};True E_{#nu}",
        lee_nbins, lee_bins, lee_nbins, lee_bins);

    TH2D* h_response_E_LEE_pi0 = new TH2D("h_response_E_LEE_pi0",
        "LEE Response Matrix;Reco E_{#nu}^{QE};True E_{#nu}",
        lee_nbins, lee_bins, lee_nbins, lee_bins);

    TH1D* h_other = new TH1D("h_other", "Other", lee_nbins, lee_bins);

    long long totalEntries = 0;

    // -------------------------------------------------------------------------
    // File loop
    // -------------------------------------------------------------------------
    for (int fileIndex = 1; fileIndex <= 10; fileIndex++)
    {
        std::stringstream ss;
        ss << OSC_MC_PREFIX << fileIndex << ".root";

        TFile* f = TFile::Open(ss.str().c_str());
        if (!f || f->IsZombie()) {
            std::cerr << "Failed to open file: " << ss.str() << std::endl;
            continue;
        }

        TTree* t = (TTree*)f->Get("MiniBooNE_CCQE");
        if (!t) {
            std::cerr << "Could not find tree MiniBooNE_CCQE in file: "
                      << ss.str() << std::endl;
            f->Close();
            continue;
        }

        int   NFSP, NUANCEChan, NuType, NuParentID;
        float Energy, RecoEnuQE, Weight, NuMomT;
        bool  PassOsc;

        std::vector<int>   *FSPType = nullptr;
        std::vector<float> *Vx = nullptr, *Vy = nullptr, *Vz = nullptr;
        std::vector<float> *MomX = nullptr, *MomY = nullptr, *MomZ = nullptr, *MomT = nullptr;

        // Branches
        t->SetBranchAddress("NFSP",      &NFSP);
        t->SetBranchAddress("FSPType",   &FSPType);
        t->SetBranchAddress("VertexX",   &Vx);
        t->SetBranchAddress("VertexY",   &Vy);
        t->SetBranchAddress("VertexZ",   &Vz);
        t->SetBranchAddress("MomX",      &MomX);
        t->SetBranchAddress("MomY",      &MomY);
        t->SetBranchAddress("MomZ",      &MomZ);
        t->SetBranchAddress("MomT",      &MomT);
        t->SetBranchAddress("NUANCEChan",&NUANCEChan);
        t->SetBranchAddress("NuType",    &NuType);
        t->SetBranchAddress("NuMomT",    &NuMomT);
        t->SetBranchAddress("NuParentID",&NuParentID);
        t->SetBranchAddress("Energy",    &Energy);
        t->SetBranchAddress("RecoEnuQE", &RecoEnuQE);
        t->SetBranchAddress("Weight",    &Weight);
        t->SetBranchAddress("PassOsc",   &PassOsc);

        int N = t->GetEntries();
        totalEntries += N;

        for (int i = 0; i < N; i++)
        {
            t->GetEntry(i);

            if (!CheckPointers(FSPType, Vx, Vy, Vz, MomX, MomY, MomZ, MomT))
                continue;

            unsigned isPi0 = sp::Pi0Details(NFSP, *FSPType, *Vx, *Vy, *Vz,
                                            *MomX, *MomY, *MomZ, *MomT);

            auto bkg_type = sp::StackHistoBkgd(false, isPi0,
                                               (sp::NuanceType_t)NUANCEChan,
                                               (sp::NuType_t)NuType,
                                               (sp::GEANT3Type_t)NuParentID);

            // Note: NuMomT, RecoEnuQE assumed in GeV already.

            switch (bkg_type) {

            case sp::kBKGD_PI0: {
            const NuanceType_t evwt = (NuanceType_t)NUANCEChan;
                // Keep only neutrino NC 1pi0 channels
                if (!(evwt == kNC1p1pi0 || evwt == kNC1n1pi0)) continue;
                h_total_true_pi0->Fill(NuMomT, Weight);
                h_total_true_LEE_pi0->Fill(NuMomT, Weight);
                if (PassOsc) {
                    h_pass_true_pi0->Fill(NuMomT, Weight);
                    h_pass_reco_pi0->Fill(RecoEnuQE, Weight);
                    h_response_pi0->Fill(RecoEnuQE, NuMomT, Weight);
                    h_pass_true_LEE_pi0->Fill(NuMomT, Weight);
                    h_pass_reco_LEE_pi0->Fill(RecoEnuQE, Weight);
                    h_response_E_LEE_pi0->Fill(RecoEnuQE, NuMomT, Weight);
                }
                break;
            }


            case sp::kBKGD_DELTA:
                if((NuanceType_t)NUANCEChan == kCC1pNg) continue;
                    h_total_true_ncdelta->Fill(NuMomT, Weight);
                    h_total_true_LEE_ncdelta->Fill(NuMomT, Weight);
                    if (PassOsc) {
                        h_pass_true_ncdelta->Fill(NuMomT, Weight);
                        h_pass_reco_ncdelta->Fill(RecoEnuQE, Weight);
                        h_response_E_ncdelta->Fill(RecoEnuQE, NuMomT, Weight);
                        // LEE binned
                        h_pass_true_LEE_ncdelta->Fill(NuMomT, Weight);
                        h_pass_reco_LEE_ncdelta->Fill(RecoEnuQE, Weight);
                        h_response_E_LEE_ncdelta->Fill(RecoEnuQE, NuMomT, Weight);
                    }
                break;
            case sp::kBKGD_NUEPIP:
                if (NuParentID == 5 && NuType == 3) {
                    h_total_true_nue->Fill(NuMomT, Weight);
                    h_total_true_LEE_nue->Fill(NuMomT, Weight);
                    if (PassOsc) {
                        h_pass_true_nue->Fill(NuMomT, Weight);
                        h_pass_reco_nue->Fill(RecoEnuQE, Weight);
                        h_response_E_nue->Fill(RecoEnuQE, NuMomT, Weight);
                        // LEE binned
                        h_pass_true_LEE_nue->Fill(NuMomT, Weight);
                        h_pass_reco_LEE_nue->Fill(RecoEnuQE, Weight);
                        h_response_E_LEE_nue->Fill(RecoEnuQE, NuMomT, Weight);
                    }
                }
                break;

            case sp::kBKGD_NUEKP:
                if (NuType == 3) {
                    h_total_true_nue->Fill(NuMomT, Weight);
                    h_total_true_LEE_nue->Fill(NuMomT, Weight);
                    if (PassOsc) {
                        h_pass_true_nue->Fill(NuMomT, Weight);
                        h_pass_reco_nue->Fill(RecoEnuQE, Weight);
                        h_response_E_nue->Fill(RecoEnuQE, NuMomT, Weight);
                        // LEE binned
                        h_pass_true_LEE_nue->Fill(NuMomT, Weight);
                        h_pass_reco_LEE_nue->Fill(RecoEnuQE, Weight);
                        h_response_E_LEE_nue->Fill(RecoEnuQE, NuMomT, Weight);
                    }
                }
                break;

            case sp::kBKGD_NUEK0:
                if (NuType == 3) {
                    h_total_true_nue->Fill(NuMomT, Weight);
                    h_total_true_LEE_nue->Fill(NuMomT, Weight);
                    if (PassOsc) {
                        h_pass_true_nue->Fill(NuMomT, Weight);
                        h_pass_reco_nue->Fill(RecoEnuQE, Weight);
                        h_response_E_nue->Fill(RecoEnuQE, NuMomT, Weight);
                        // LEE binned
                        h_pass_true_LEE_nue->Fill(NuMomT, Weight);
                        h_pass_reco_LEE_nue->Fill(RecoEnuQE, Weight);
                        h_response_E_LEE_nue->Fill(RecoEnuQE, NuMomT, Weight);
                    }
                }
                break;

            default:
                h_other->Fill(RecoEnuQE, Weight);
                break;
            }
        } // event loop

        f->Close();
    } // file loop

    // ========================================================================
    // Clone raw (pre-normalisation) smearing matrices
    // ========================================================================
    TH2D* h_smearing_E_nue         = (TH2D*)h_response_E_nue->Clone("h_smearing_E_nue");
    TH2D* h_smearing_E_ncdelta     = (TH2D*)h_response_E_ncdelta->Clone("h_smearing_E_ncdelta");
    TH2D* h_smearing_pi0           = (TH2D*)h_response_pi0->Clone("h_smearing_pi0");

    TH2D* h_smearing_E_LEE_nue     = (TH2D*)h_response_E_LEE_nue->Clone("h_smearing_E_LEE_nue");
    TH2D* h_smearing_E_LEE_ncdelta = (TH2D*)h_response_E_LEE_ncdelta->Clone("h_smearing_E_LEE_ncdelta");
    TH2D* h_smearing_E_LEE_pi0     = (TH2D*)h_response_E_LEE_pi0->Clone("h_smearing_E_LEE_pi0");

    // -------------------------------------------------------------------------
    // Normalise response matrices (per true bin)
    // -------------------------------------------------------------------------
    for (int iy = 1; iy <= nbins_true_nue; iy++) {
        double den = h_total_true_nue->GetBinContent(iy);
        if (den <= 0) continue;
        for (int ix = 1; ix <= nbins_reco_nue; ix++) {
            double val = h_response_E_nue->GetBinContent(ix, iy) / den;
            h_response_E_nue->SetBinContent(ix, iy, val);
        }
    }

    for (int iy = 1; iy <= nbins_true_ncdelta; iy++) {
        double den = h_total_true_ncdelta->GetBinContent(iy);
        if (den <= 0) continue;
        for (int ix = 1; ix <= nbins_reco_ncdelta; ix++) {
            double val = h_response_E_ncdelta->GetBinContent(ix, iy) / den;
            h_response_E_ncdelta->SetBinContent(ix, iy, val);
        }
    }

    for (int iy = 1; iy <= nbins_true_pi0; iy++) {
        double den = h_total_true_pi0->GetBinContent(iy);
        if (den <= 0) continue;
        for (int ix = 1; ix <= nbins_reco_pi0; ix++) {
            double val = h_response_pi0->GetBinContent(ix, iy) / den;
            h_response_pi0->SetBinContent(ix, iy, val);
        }
    }

    // LEE-binned normalisation
    auto NormLEE = [&](TH1D* htot, TH2D* hresp){
        for (int iy = 1; iy <= lee_nbins; iy++) {
            double den = htot->GetBinContent(iy);
            if (den <= 0) continue;
            for (int ix = 1; ix <= lee_nbins; ix++) {
                double val = hresp->GetBinContent(ix, iy) / den;
                hresp->SetBinContent(ix, iy, val);
            }
        }
    };

    NormLEE(h_total_true_LEE_nue,     h_response_E_LEE_nue);
    NormLEE(h_total_true_LEE_ncdelta, h_response_E_LEE_ncdelta);
    NormLEE(h_total_true_LEE_pi0,     h_response_E_LEE_pi0);

    // -------------------------------------------------------------------------
    // Indexed response matrices (normalised)
    // -------------------------------------------------------------------------
    TH2D* h_response_E_idx_nue =
        new TH2D("h_response_E_idx_nue",";Reco E bin;True E bin",
                 nbins_reco_nue, 0, nbins_reco_nue,
                 nbins_true_nue, 0, nbins_true_nue);
    CopyAndLabelIndexed(h_response_E_nue, h_response_E_idx_nue);

    TH2D* h_response_E_idx_ncdelta =
        new TH2D("h_response_E_idx_ncdelta",";Reco E bin;True E bin",
                 nbins_reco_ncdelta, 0, nbins_reco_ncdelta,
                 nbins_true_ncdelta, 0, nbins_true_ncdelta);
    CopyAndLabelIndexed(h_response_E_ncdelta, h_response_E_idx_ncdelta);

    TH2D* h_response_E_idx_pi0 =
        new TH2D("h_response_E_idx_pi0",";Reco E bin;True E bin",
                 nbins_reco_pi0, 0, nbins_reco_pi0,
                 nbins_true_pi0, 0, nbins_true_pi0);
    CopyAndLabelIndexed(h_response_pi0, h_response_E_idx_pi0);

    // -------------------------------------------------------------------------
    // Indexed RAW smearing matrices (standard binning)
    // -------------------------------------------------------------------------
    TH2D* h_smearing_E_idx_nue =
        new TH2D("smearing_nue_idx",";Reco E bin;True E bin",
                 nbins_reco_nue, 0, nbins_reco_nue,
                 nbins_true_nue, 0, nbins_true_nue);
    CopyAndLabelIndexed(h_smearing_E_nue, h_smearing_E_idx_nue);

    TH2D* h_smearing_E_idx_ncdelta =
        new TH2D("smearing_ncdelta_idx",";Reco E bin;True E bin",
                 nbins_reco_ncdelta, 0, nbins_reco_ncdelta,
                 nbins_true_ncdelta, 0, nbins_true_ncdelta);
    CopyAndLabelIndexed(h_smearing_E_ncdelta, h_smearing_E_idx_ncdelta);

    TH2D* h_smearing_E_idx_pi0 =
        new TH2D("smearing_pi0_idx",";Reco E bin;True E bin",
                 nbins_reco_pi0, 0, nbins_reco_pi0,
                 nbins_true_pi0, 0, nbins_true_pi0);
    CopyAndLabelIndexed(h_smearing_pi0, h_smearing_E_idx_pi0);
    
    // -------------------------------------------------------------------------
    // Events per MeV histograms
    // -------------------------------------------------------------------------
    auto h_total_true_nue_EpMeV     = (TH1D*)h_total_true_nue->Clone("h_total_true_nue_EpMeV");
    auto h_pass_true_nue_EpMeV      = (TH1D*)h_pass_true_nue->Clone("h_pass_true_nue_EpMeV");
    auto h_pass_reco_nue_EpMeV      = (TH1D*)h_pass_reco_nue->Clone("h_pass_reco_nue_EpMeV");

    auto h_total_true_ncdelta_EpMeV = (TH1D*)h_total_true_ncdelta->Clone("h_total_true_ncdelta_EpMeV");
    auto h_pass_true_ncdelta_EpMeV  = (TH1D*)h_pass_true_ncdelta->Clone("h_pass_true_ncdelta_EpMeV");
    auto h_pass_reco_ncdelta_EpMeV  = (TH1D*)h_pass_reco_ncdelta->Clone("h_pass_reco_ncdelta_EpMeV");

    auto h_total_true_pi0_EpMeV     = (TH1D*)h_total_true_pi0->Clone("h_total_true_pi0_EpMeV");
    auto h_pass_true_pi0_EpMeV      = (TH1D*)h_pass_true_pi0->Clone("h_pass_true_pi0_EpMeV");
    auto h_pass_reco_pi0_EpMeV      = (TH1D*)h_pass_reco_pi0->Clone("h_pass_reco_pi0_EpMeV");

    NormaliseByBinWidth(h_total_true_nue_EpMeV);
    NormaliseByBinWidth(h_pass_true_nue_EpMeV);
    NormaliseByBinWidth(h_pass_reco_nue_EpMeV);

    NormaliseByBinWidth(h_total_true_ncdelta_EpMeV);
    NormaliseByBinWidth(h_pass_true_ncdelta_EpMeV);
    NormaliseByBinWidth(h_pass_reco_ncdelta_EpMeV);

    NormaliseByBinWidth(h_total_true_pi0_EpMeV);
    NormaliseByBinWidth(h_pass_true_pi0_EpMeV);
    NormaliseByBinWidth(h_pass_reco_pi0_EpMeV);

    // -------------------------------------------------------------------------
    // Output directories
    // -------------------------------------------------------------------------
   const char* png_dir  = MATRIX_PNG_OUTPUT_DIR;
   const char* root_dir = MATRIX_ROOT_OUTPUT_DIR;

    // -------------------------------------------------------------------------
    // Plotting helpers
    // -------------------------------------------------------------------------
    auto SaveMatrix = [&](TH2* h, const char* fname){
        TCanvas* c = new TCanvas(fname, fname, 900, 800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.08);
        h->LabelsOption("h");
        h->Draw("COLZ");
        std::string out = std::string(png_dir) + "/" + fname + "_all.png";

        c->SaveAs(out.c_str());
        delete c;
    };

    auto Save1D = [&](TH1* h, const char* fname){
        TCanvas* c = new TCanvas(fname, fname, 800, 800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.15);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.08);
        h->Draw("HIST E");
        std::string out = std::string(png_dir) + "/" + fname + "_all.png";
        c->SaveAs(out.c_str());
        delete c;
    };

    // -------------------------------------------------------------------------
    // Save PNGs
    // -------------------------------------------------------------------------
    // Normalised response matrices
    SaveMatrix(h_response_E_nue,        "Nue_ResponseMatrix_Energy");
    SaveMatrix(h_response_E_idx_nue,    "Nue_ResponseMatrix_Energy_Indexed");
    SaveMatrix(h_response_E_LEE_nue,    "Nue_ResponseMatrix_Energy_LEE");

    SaveMatrix(h_response_E_ncdelta,    "NCDelta_ResponseMatrix_Energy");
    SaveMatrix(h_response_E_idx_ncdelta,"NCDelta_ResponseMatrix_Energy_Indexed");
    SaveMatrix(h_response_E_LEE_ncdelta,"NCDelta_ResponseMatrix_Energy_LEE");

    SaveMatrix(h_response_pi0,          "Pi0_ResponseMatrix_Energy");
    SaveMatrix(h_response_E_idx_pi0,    "Pi0_ResponseMatrix_Energy_Indexed");
    SaveMatrix(h_response_E_LEE_pi0,    "Pi0_ResponseMatrix_Energy_LEE");

    // Indexed RAW matrices (standard binning)
    SaveMatrix(h_smearing_E_idx_nue,     "Nue_SmearingMatrix_Energy_Indexed");
    SaveMatrix(h_smearing_E_idx_ncdelta, "NCDelta_SmearingMatrix_Energy_Indexed");
    SaveMatrix(h_smearing_E_idx_pi0,     "Pi0_SmearingMatrix_Energy_Indexed");

    // 1D spectra
    Save1D(h_total_true_nue,     "Nue_TrueEnergy_Total");
    Save1D(h_pass_true_nue,      "Nue_TrueEnergy_Passed");
    Save1D(h_pass_reco_nue,      "Nue_RecoEnergy_Passed");

    Save1D(h_total_true_ncdelta, "NCDelta_TrueEnergy_Total");
    Save1D(h_pass_true_ncdelta,  "NCDelta_TrueEnergy_Passed");
    Save1D(h_pass_reco_ncdelta,  "NCDelta_RecoEnergy_Passed");

    Save1D(h_total_true_pi0,     "Pi0_TrueEnergy_Total");
    Save1D(h_pass_true_pi0,      "Pi0_TrueEnergy_Passed");
    Save1D(h_pass_reco_pi0,      "Pi0_RecoEnergy_Passed");

    // EpMeV spectra
    Save1D(h_total_true_nue_EpMeV,     "Nue_TrueEnergy_Total_EpMeV");
    Save1D(h_pass_true_nue_EpMeV,      "Nue_TrueEnergy_Passed_EpMeV");
    Save1D(h_pass_reco_nue_EpMeV,      "Nue_RecoEnergy_Passed_EpMeV");

    Save1D(h_total_true_ncdelta_EpMeV, "NCDelta_TrueEnergy_Total_EpMeV");
    Save1D(h_pass_true_ncdelta_EpMeV,  "NCDelta_TrueEnergy_Passed_EpMeV");
    Save1D(h_pass_reco_ncdelta_EpMeV,  "NCDelta_RecoEnergy_Passed_EpMeV");

    Save1D(h_total_true_pi0_EpMeV,     "Pi0_TrueEnergy_Total_EpMeV");
    Save1D(h_pass_true_pi0_EpMeV,      "Pi0_TrueEnergy_Passed_EpMeV");
    Save1D(h_pass_reco_pi0_EpMeV,      "Pi0_RecoEnergy_Passed_EpMeV");

    // -------------------------------------------------------------------------
    // ROOT output files
    // -------------------------------------------------------------------------
    {
        std::string outfile1 = std::string(root_dir) + "/all.root";
        TFile out1(outfile1.c_str(), "RECREATE");

        h_total_true_nue->Write("total_true_nue");
        h_total_true_ncdelta->Write("total_true_ncdelta");
        h_total_true_pi0->Write("total_true_pi0");

        h_pass_true_nue->Write("pass_true_nue");
        h_pass_true_ncdelta->Write("pass_true_ncdelta");
        h_pass_true_pi0->Write("pass_true_pi0");

        h_pass_reco_nue->Write("pass_reco_nue");
        h_pass_reco_ncdelta->Write("pass_reco_ncdelta");
        h_pass_reco_pi0->Write("pass_reco_pi0");

        // Normalised response matrices
        h_response_E_nue->Write("response_nue");
        h_response_E_ncdelta->Write("response_ncdelta");
        h_response_pi0->Write("response_pi0");

        // Raw (pre-normalisation) response matrices
        h_smearing_E_nue->Write("smearing_nue");
        h_smearing_E_ncdelta->Write("smearing_ncdelta");
        h_smearing_pi0->Write("smearing_pi0");

        // Indexed RAW (standard-binning) matrices
        h_smearing_E_idx_nue->Write("smearing_nue_idx");
        h_smearing_E_idx_ncdelta->Write("smearing_ncdelta_idx");
        h_smearing_E_idx_pi0->Write("smearing_pi0_idx");

        out1.Close();

        std::string outfile2 = std::string(root_dir) + "/response_LEE_all.root";
        TFile out2(outfile2.c_str(), "RECREATE");

        h_total_true_LEE_nue->Write("total_true_LEE_nue");
        h_total_true_LEE_ncdelta->Write("total_true_LEE_ncdelta");
        h_total_true_LEE_pi0->Write("total_true_LEE_pi0");

        h_pass_true_LEE_nue->Write("pass_true_LEE_nue");
        h_pass_true_LEE_ncdelta->Write("pass_true_LEE_ncdelta");
        h_pass_true_LEE_pi0->Write("pass_true_LEE_pi0");

        h_pass_reco_LEE_nue->Write("pass_reco_LEE_nue");
        h_pass_reco_LEE_ncdelta->Write("pass_reco_LEE_ncdelta");
        h_pass_reco_LEE_pi0->Write("pass_reco_LEE_pi0");

        h_response_E_LEE_nue->Write("response_LEE_nue");
        h_response_E_LEE_ncdelta->Write("response_LEE_ncdelta");
        h_response_E_LEE_pi0->Write("response_LEE_pi0");

        // Raw LEE response matrices (no indexed LEE raw requested)
        h_smearing_E_LEE_nue->Write("smearing_LEE_nue");
        h_smearing_E_LEE_ncdelta->Write("smearing_LEE_ncdelta");
        h_smearing_E_LEE_pi0->Write("smearing_LEE_pi0");

        out2.Close();
    }

    std::cout << "Done!" << std::endl;
}