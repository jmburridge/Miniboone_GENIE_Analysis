// Pi0ResponseMatrix_indexed.C
// Builds π0 response matrices (energy & momentum) and classifies π0 events
// Now also produces a separate MiniBooNE LEE-binned energy response matrix
// in a separate ROOT file.
//'Technote-style' binning scheme originally from MiniBooNE technote 214: https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=9914&filename=unfolding_v2.0.pdf&version=2
//MinibooNE LEE binning scheme from MiniBooNE LEE analysis: https://arxiv.org/abs/1805.12028

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"

#include "MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedTypes.h"
#include "MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.h"
#include "MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.cxx"

using namespace sp;

// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------
bool CheckPointers(std::vector<int>* FSPType, std::vector<float>* Vx, std::vector<float>* Vy,
                   std::vector<float>* Vz, std::vector<float>* MomX, std::vector<float>* MomY,
                   std::vector<float>* MomZ, std::vector<float>* MomT)
{
    return (FSPType && Vx && Vy && Vz && MomX && MomY && MomZ && MomT);
}

void CopyAndLabelIndexed(TH2D* source, TH2D* target)
{
    int nx = source->GetNbinsX();
    int ny = source->GetNbinsY();

    // Copy contents
    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            target->SetBinContent(ix, iy, source->GetBinContent(ix, iy));
        }
    }

    // Label axes with bin indices
    for (int i = 1; i <= nx; ++i)
        target->GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
    for (int j = 1; j <= ny; ++j)
        target->GetYaxis()->SetBinLabel(j, std::to_string(j).c_str());
}

// -----------------------------------------------------------------------------
// Main macro
// -----------------------------------------------------------------------------
void Pi0ResponseMatrix()
{
    // -------------------------------------------------------------------------
    // Binning schemes
    // -------------------------------------------------------------------------
    // Original 9-bin scheme
    const int nbins = 9;
    double recoE_bins[nbins + 1] = {0.150, 0.200, 0.300, 0.375, 0.475, 0.550, 0.675, 0.800, 1.000, 1.200};
    double trueE_bins[nbins + 1] = {0.250, 0.500, 0.600, 0.700, 0.900, 1.100, 1.250, 1.500, 2.00, 3.00};
    double trueP_bins[nbins + 1] = {0.050, 0.200, 0.300, 0.400, 0.500, 0.700, 0.750, 0.800, 1.00, 1.200};
    double recoP_bins[nbins + 1] = {0.100, 0.300, 0.450, 0.500, 0.750, 0.800, 0.900, 1.000, 1.100, 1.200};

    // MiniBooNE LEE 11-bin scheme (square matrix: true and reco use same bins)
    const int lee_nbins = 11;
    double lee_bins[lee_nbins + 1] = {0.2, 0.3, 0.375, 0.475, 0.55, 0.675, 0.8, 0.95, 1.1, 1.3, 1.5, 3.0};

    // -------------------------------------------------------------------------
    // Histograms: Technote-style 9-bin scheme: energy & momentum response matrices
    // -------------------------------------------------------------------------
    TH1D *h_total_true       = new TH1D("h_total_true",       "Total True E_{#nu};E_{#nu} [GeV];Events",               nbins, trueE_bins);
    TH1D *h_total_true_mom   = new TH1D("h_total_true_mom",   "Total True Pi0 Momentum;p_{Pi0} [GeV/c];Events",        nbins, trueP_bins);
    TH1D *h_pass_true        = new TH1D("h_pass_true",        "Passed Events by True E_{#nu};E_{#nu} [GeV];Events",    nbins, trueE_bins);
    TH1D *h_pass_true_mom    = new TH1D("h_pass_true_mom",    "Passed Events by True p_{Pi0};p_{Pi0} [GeV/c];Events",  nbins, trueP_bins);
    TH1D *h_pass_reco_energy = new TH1D("h_pass_reco_energy", "Passed Events by Reco E_{#nu}^{QE};Reco E_{#nu}^{QE} [GeV];Events",
                                        nbins, recoE_bins);

    TH2D *h_response_E = new TH2D("h_response_E",
                                  "Reco vs True E_{#nu};Reco E_{#nu}^{QE} [GeV];True E_{#nu} [GeV]",
                                  nbins, recoE_bins, nbins, trueE_bins);
    TH2D *h_response_P = new TH2D("h_response_P",
                                  "Reco vs True p_{Pi0};Reco E_{#nu}^{QE} [GeV];True p_{Pi0} [GeV/c]",
                                  nbins, recoP_bins, nbins, trueP_bins);

    TH2D *h_response_E_CC = new TH2D("h_response_E_CC",
                                     "Reco vs True E_{#nu} (CC);Reco E_{#nu}^{QE} [GeV];True E_{#nu} [GeV]",
                                     nbins, recoE_bins, nbins, trueE_bins);
    TH2D *h_response_E_NC = new TH2D("h_response_E_NC",
                                     "Reco vs True E_{#nu} (NC);Reco E_{#nu}^{QE} [GeV];True E_{#nu} [GeV]",
                                     nbins, recoE_bins, nbins, trueE_bins);
    TH2D *h_response_P_CC = new TH2D("h_response_P_CC",
                                     "Reco vs True p_{#pi^{0}} (CC);Reco E_{#nu}^{QE} [GeV];True p_{#pi^{0}} [GeV/c]",
                                     nbins, recoP_bins, nbins, trueP_bins);
    TH2D *h_response_P_NC = new TH2D("h_response_P_NC",
                                     "Reco vs True p_{#pi^{0}} (NC);Reco E_{#nu}^{QE} [GeV];True p_{#pi^{0}} [GeV/c]",
                                     nbins, recoP_bins, nbins, trueP_bins);

    TH1D *h_total_true_mom_CC = new TH1D("h_total_true_mom_CC",
                                         "Total True Pi0 Momentum (CC);p_{#pi^{0}} [GeV/c];Events",
                                         nbins, trueP_bins);
    TH1D *h_total_true_mom_NC = new TH1D("h_total_true_mom_NC",
                                         "Total True Pi0 Momentum (NC);p_{#pi^{0}} [GeV/c];Events",
                                         nbins, trueP_bins);
    TH1D *h_eff_pi0_mom_CC    = new TH1D("h_eff_pi0_mom_CC",
                                         "Pi0 Efficiency (CC);True p_{#pi^{0}} [GeV/c];Efficiency",
                                         nbins, trueP_bins);
    TH1D *h_eff_pi0_mom_NC    = new TH1D("h_eff_pi0_mom_NC",
                                         "Pi0 Efficiency (NC);True p_{#pi^{0}} [GeV/c];Efficiency",
                                         nbins, trueP_bins);

    TH1D *h_eff_pi0_mom = new TH1D("h_eff_pi0_mom",
                                   "Pi0 Efficiency;True Pi0} [GeV/c];Efficiency",
                                   nbins, trueP_bins);
    TH1D *h_eff_pi0     = new TH1D("h_eff_pi0",
                                   " ;True #pi_{0} Energy [GeV];Efficiency",
                                   nbins, trueE_bins);

    TH1D *h_ratio_passRecoE_totalTrueP =
        new TH1D("h_ratio_passRecoE_totalTrueP",
                 "Ratio: Passed Reco E^{QE} / Total True p_{#pi^{0}};Bin;Ratio",
                 nbins, 1, nbins + 1);

    TH1F *h_pi0_nuance_class =
        new TH1F("h_pi0_nuance_class", "Pi0 by NUANCEChan;Interaction Type;Events", 7, 0.5, 7.5);

    const char* nuance_labels[] = {
        "CC res (4)",
        "NC res (6)",
        "NC res (8)",
        "CC #Delta^{++} (18)",
        "CC #Delta^{+} (19)",
        "NC #Delta^{+} (22)",
        "NC #Delta^{0} (26)"
    };
    for (int i = 0; i < 7; ++i)
        h_pi0_nuance_class->GetXaxis()->SetBinLabel(i + 1, nuance_labels[i]);

    // Histogram of event counts: CC reco, NC reco, CC total, CC passed, NC total, NC passed
    TH1F* h_pi0_event_counts =
        new TH1F("h_pi0_event_counts", "Pi0 Event Counts;Category;Events", 6, 0.5, 6.5);
    h_pi0_event_counts->GetXaxis()->SetBinLabel(1, "CC Reco");
    h_pi0_event_counts->GetXaxis()->SetBinLabel(2, "NC Reco");
    h_pi0_event_counts->GetXaxis()->SetBinLabel(3, "CC Total");
    h_pi0_event_counts->GetXaxis()->SetBinLabel(4, "CC Passed");
    h_pi0_event_counts->GetXaxis()->SetBinLabel(5, "NC Total");
    h_pi0_event_counts->GetXaxis()->SetBinLabel(6, "NC Passed");

    // -------------------------------------------------------------------------
    // MiniBooNE LEE-binned energy histograms (square matrix)
    // -------------------------------------------------------------------------
    TH1D *h_total_true_LEE = new TH1D("h_total_true_LEE",
                                      "Total True E_{#nu} (MiniBooNE LEE);E_{#nu} [GeV];Events",
                                      lee_nbins, lee_bins);
    TH1D *h_pass_true_LEE  = new TH1D("h_pass_true_LEE",
                                      "Passed Events by True E_{#nu} (MiniBooNE LEE);E_{#nu} [GeV];Events",
                                      lee_nbins, lee_bins);

    TH2D *h_response_E_LEE = new TH2D("h_response_E_LEE",
                                      "Reco vs True E_{#nu} (MiniBooNE LEE);Reco E_{#nu}^{QE} [GeV];True E_{#nu} [GeV]",
                                      lee_nbins, lee_bins,
                                      lee_nbins, lee_bins);

    // -------------------------------------------------------------------------
    // Counters
    // -------------------------------------------------------------------------
    long long totalEntries   = 0;
    long long pi0count       = 0;
    long long count_total_CC = 0;
    long long count_total_NC = 0;
    long long count_passed_CC = 0;
    long long count_passed_NC = 0;
    long long count_reco_CC   = 0;
    long long count_reco_NC   = 0;

    // -------------------------------------------------------------------------
    // Loop over input ROOT files
    // -------------------------------------------------------------------------
    for (int fileIndex = 1; fileIndex <= 10; ++fileIndex) {

        std::stringstream filename;
        filename << "MiniBooNEDatasets2023/output_osc_mc_detail_" << fileIndex << ".root";
        TFile* file = TFile::Open(filename.str().c_str());
        if (!file || file->IsZombie()) {
            continue;
        }

        TTree* tree = (TTree*)file->Get("MiniBooNE_CCQE");
        if (!tree) {
            file->Close();
            continue;
        }

        int   NFSP, NUANCEChan, NuType, NuParentID;
        float Energy, RecoEnuQE, Weight, NuMomT;
        bool  PassOsc;

        std::vector<int>   *FSPType = nullptr;
        std::vector<float> *Vx = nullptr, *Vy = nullptr, *Vz = nullptr;
        std::vector<float> *MomX = nullptr, *MomY = nullptr, *MomZ = nullptr, *MomT = nullptr;

        tree->SetBranchAddress("NFSP",       &NFSP);
        tree->SetBranchAddress("FSPType",    &FSPType);
        tree->SetBranchAddress("VertexX",    &Vx);
        tree->SetBranchAddress("VertexY",    &Vy);
        tree->SetBranchAddress("VertexZ",    &Vz);
        tree->SetBranchAddress("MomX",       &MomX);
        tree->SetBranchAddress("MomY",       &MomY);
        tree->SetBranchAddress("MomZ",       &MomZ);
        tree->SetBranchAddress("MomT",       &MomT);
        tree->SetBranchAddress("NUANCEChan", &NUANCEChan);
        tree->SetBranchAddress("NuType",     &NuType);
        tree->SetBranchAddress("NuMomT",     &NuMomT);
        tree->SetBranchAddress("NuParentID", &NuParentID);
        tree->SetBranchAddress("Energy",     &Energy);
        tree->SetBranchAddress("RecoEnuQE",  &RecoEnuQE);
        tree->SetBranchAddress("Weight",     &Weight);
        tree->SetBranchAddress("PassOsc",    &PassOsc);

        int nentries = tree->GetEntries();
        totalEntries += nentries;

        // ---------------------------------------------------------------------
        // Event loop
        // ---------------------------------------------------------------------
        for (int i = 0; i < nentries; ++i) {

            tree->GetEntry(i);
            if (!CheckPointers(FSPType, Vx, Vy, Vz, MomX, MomY, MomZ, MomT)) {
                continue;
            }

            unsigned npi0 = sp::Pi0Details(NFSP, *FSPType, *Vx, *Vy, *Vz,
                                           *MomX, *MomY, *MomZ, *MomT);
            bool Event_is_pi0 = (npi0 > 0);

            StackedBkgdType_t bkgd_type = StackHistoBkgd(
                false,                        // Event_is_dirt (FIXME in future)
                Event_is_pi0,                 // Event_is_pi0
                static_cast<NuanceType_t>(NUANCEChan),
                static_cast<NuType_t>(NuType),
                static_cast<GEANT3Type_t>(NuParentID)
            );

            if (bkgd_type != kBKGD_PI0) {
                continue;
            }

            // Find π0 for momentum
            int pi0_idx = -1;
            for (size_t j = 0; j < FSPType->size(); ++j) {
                if (FSPType->at(j) == kPION0) {
                    pi0_idx = (int)j;
                    break;
                }
            }
            if (pi0_idx == -1) {
                continue;
            }

            // NUANCE classification
            switch (NUANCEChan) {
                case 4:  h_pi0_nuance_class->Fill(1); break;
                case 6:  h_pi0_nuance_class->Fill(2); break;
                case 8:  h_pi0_nuance_class->Fill(3); break;
                case 18: h_pi0_nuance_class->Fill(4); break;
                case 19: h_pi0_nuance_class->Fill(5); break;
                case 22: h_pi0_nuance_class->Fill(6); break;
                case 26: h_pi0_nuance_class->Fill(7); break;
                default: continue;
            }

            ++pi0count;

            double p_pi0 = MomT->at(pi0_idx);  // using MomT as in original code
            bool isCC = (NUANCEChan == 4  || NUANCEChan == 18 || NUANCEChan == 19);
            bool isNC = (NUANCEChan == 6  || NUANCEChan == 8  ||
                         NUANCEChan == 22 || NUANCEChan == 26);

            double pot_weight = 1.0;  // placeholder scaling, same as original

            // Fill truth distributions (original binning)
            h_total_true->Fill(NuMomT, Weight);
            h_total_true_mom->Fill(p_pi0, Weight);

            if (isCC) {
                h_total_true_mom_CC->Fill(p_pi0, Weight * pot_weight);
                ++count_total_CC;
            }
            if (isNC) {
                h_total_true_mom_NC->Fill(p_pi0, Weight * pot_weight);
                ++count_total_NC;
            }

            // Fill truth distributions with LEE binning
            h_total_true_LEE->Fill(NuMomT, Weight);

            if (PassOsc) {
                // Original binning
                h_pass_true->Fill(NuMomT, Weight * pot_weight);
                h_pass_true_mom->Fill(p_pi0, Weight * pot_weight);
                h_response_E->Fill(RecoEnuQE, NuMomT, Weight * pot_weight);
                h_response_P->Fill(RecoEnuQE, p_pi0,  Weight * pot_weight);
                h_pass_reco_energy->Fill(RecoEnuQE, Weight * pot_weight);

                // LEE-binned histograms
                h_pass_true_LEE->Fill(NuMomT, Weight * pot_weight);
                h_response_E_LEE->Fill(RecoEnuQE, NuMomT, Weight * pot_weight);

                if (isCC) {
                    h_response_P_CC->Fill(RecoEnuQE, p_pi0, Weight * pot_weight);
                    h_response_E_CC->Fill(RecoEnuQE, NuMomT, Weight * pot_weight);
                    h_eff_pi0_mom_CC->Fill(p_pi0, Weight * pot_weight);
                    ++count_passed_CC;
                    ++count_reco_CC;
                } else if (isNC) {
                    h_response_P_NC->Fill(RecoEnuQE, p_pi0, Weight * pot_weight);
                    h_response_E_NC->Fill(RecoEnuQE, NuMomT, Weight * pot_weight);
                    h_eff_pi0_mom_NC->Fill(p_pi0, Weight * pot_weight);
                    ++count_passed_NC;
                    ++count_reco_NC;
                }
            }
        }

        file->Close();
    }

    // -------------------------------------------------------------------------
    // Normalise thechnote 9-bin response matrices by true distributions
    // -------------------------------------------------------------------------
    for (int iy = 1; iy <= nbins; ++iy) {
        double E_true = h_total_true->GetBinContent(iy);
        double P_true = h_total_true_mom->GetBinContent(iy);

        if (E_true > 0.0) {
            for (int ix = 1; ix <= nbins; ++ix) {
                double val = h_response_E->GetBinContent(ix, iy);
                h_response_E->SetBinContent(ix, iy, val / E_true);
            }
        }

        if (P_true > 0.0) {
            for (int ix = 1; ix <= nbins; ++ix) {
                double val = h_response_P->GetBinContent(ix, iy);
                h_response_P->SetBinContent(ix, iy, val / P_true);
            }
        }
    }

    // Ratio histogram: passed true E / total true p (original code)
    for (int i = 1; i <= nbins; ++i) {
        double recoE = h_pass_true->GetBinContent(i);
        double trueP = h_total_true_mom->GetBinContent(i);
        if (trueP > 0.0)
            h_ratio_passRecoE_totalTrueP->SetBinContent(i, recoE / trueP);
        else
            h_ratio_passRecoE_totalTrueP->SetBinContent(i, 0.0);
    }

    // CC and NC response matrices by true p
    for (int iy = 1; iy <= nbins; ++iy) {
        double P_CC = h_total_true_mom_CC->GetBinContent(iy);
        double P_NC = h_total_true_mom_NC->GetBinContent(iy);

        if (P_CC > 0.0) {
            for (int ix = 1; ix <= nbins; ++ix) {
                double val = h_response_P_CC->GetBinContent(ix, iy);
                h_response_P_CC->SetBinContent(ix, iy, val / P_CC);
            }
        }

        if (P_NC > 0.0) {
            for (int ix = 1; ix <= nbins; ++ix) {
                double val = h_response_P_NC->GetBinContent(ix, iy);
                h_response_P_NC->SetBinContent(ix, iy, val / P_NC);
            }
        }
    }

    // -------------------------------------------------------------------------
    // Normalise LEE-binned energy response matrix (square 11x11)
    // -------------------------------------------------------------------------
    for (int iy = 1; iy <= lee_nbins; ++iy) {
        double E_true_LEE = h_total_true_LEE->GetBinContent(iy);
        if (E_true_LEE <= 0.0) continue;

        for (int ix = 1; ix <= lee_nbins; ++ix) {
            double val = h_response_E_LEE->GetBinContent(ix, iy);
            h_response_E_LEE->SetBinContent(ix, iy, val / E_true_LEE);
        }
    }

    // -------------------------------------------------------------------------
    // Total π0 efficiency with binomial errors (original logic)
    // -------------------------------------------------------------------------
    for (int i = 1; i <= nbins; ++i) {
        double num        = h_pass_true_mom->GetBinContent(i);
        double den        = h_total_true_mom->GetBinContent(i);
        double num_energy = h_pass_true->GetBinContent(i);
        double den_energy = h_total_true->GetBinContent(i);

        if (den > 0.0) {
            double eff        = num / den;
            double err        = std::sqrt(eff * (1.0 - eff) / den);
            double eff_energy = (den_energy > 0.0) ? num_energy / den_energy : 0.0;
            // NOTE: the original code uses 'eff' and 'den' for err_energy as well
            double err_energy = std::sqrt(eff * (1.0 - eff) / den);

            h_eff_pi0_mom->SetBinContent(i, eff);
            h_eff_pi0_mom->SetBinError(i, err);

            h_eff_pi0->SetBinContent(i, eff_energy);
            h_eff_pi0->SetBinError(i, err_energy);
        } else {
            h_eff_pi0_mom->SetBinContent(i, 0.0);
            h_eff_pi0_mom->SetBinError(i, 0.0);
            h_eff_pi0->SetBinContent(i, 0.0);
            h_eff_pi0->SetBinError(i, 0.0);
        }
    }

    // CC π0 efficiency with binomial errors
    for (int i = 1; i <= nbins; ++i) {
        double num = h_eff_pi0_mom_CC->GetBinContent(i); // numerator stored here
        double den = h_total_true_mom_CC->GetBinContent(i);

        if (den > 0.0) {
            double eff = num / den;
            double err = std::sqrt(eff * (1.0 - eff) / den);
            h_eff_pi0_mom_CC->SetBinContent(i, eff);
            h_eff_pi0_mom_CC->SetBinError(i, err);
        } else {
            h_eff_pi0_mom_CC->SetBinContent(i, 0.0);
            h_eff_pi0_mom_CC->SetBinError(i, 0.0);
        }
    }

    // NC π0 efficiency with binomial errors
    for (int i = 1; i <= nbins; ++i) {
        double num = h_eff_pi0_mom_NC->GetBinContent(i);
        double den = h_total_true_mom_NC->GetBinContent(i);

        if (den > 0.0) {
            double eff = num / den;
            double err = std::sqrt(eff * (1.0 - eff) / den);
            h_eff_pi0_mom_NC->SetBinContent(i, eff);
            h_eff_pi0_mom_NC->SetBinError(i, err);
        } else {
            h_eff_pi0_mom_NC->SetBinContent(i, 0.0);
            h_eff_pi0_mom_NC->SetBinError(i, 0.0);
        }
    }

    // -------------------------------------------------------------------------
    // Fill event-count histogram from existing histograms
    // -------------------------------------------------------------------------
    h_pi0_event_counts->SetBinContent(1, h_response_P_CC->GetEntries());   // CC Reco
    h_pi0_event_counts->SetBinContent(2, h_response_P_NC->GetEntries());   // NC Reco
    h_pi0_event_counts->SetBinContent(3, h_total_true_mom_CC->Integral()); // CC Total
    h_pi0_event_counts->SetBinContent(4, h_eff_pi0_mom_CC->Integral());    // CC Passed
    h_pi0_event_counts->SetBinContent(5, h_total_true_mom_NC->Integral()); // NC Total
    h_pi0_event_counts->SetBinContent(6, h_eff_pi0_mom_NC->Integral());    // NC Passed

    // -------------------------------------------------------------------------
    // Indexed copies of original response matrices
    // -------------------------------------------------------------------------
    TH2D* h_response_E_idx = new TH2D("h_response_E_idx", ";Reco E bin;True E bin",
                                      nbins, 0, nbins, nbins, 0, nbins);
    TH2D* h_response_P_idx = new TH2D("h_response_P_idx", ";Reco E bin;True p_{#pi^{0}} bin",
                                      nbins, 0, nbins, nbins, 0, nbins);

    CopyAndLabelIndexed(h_response_E, h_response_E_idx);
    CopyAndLabelIndexed(h_response_P, h_response_P_idx);

    // -------------------------------------------------------------------------
    // ROOT style
    // -------------------------------------------------------------------------
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetOptStat(0);

    auto DrawMatrix = [&](TH2* hist, const std::string& name)
    {
        TCanvas* c = new TCanvas(name.c_str(), name.c_str(), 900, 800);
        hist->LabelsOption("h");
        hist->Draw("COLZ");
        c->SaveAs(Form("Pi0Analysis/%s.png", name.c_str()));
    };

    auto PlotEventsAndEfficiency =
        [&](TH1* h_true_pass, TH1* h_reco_pass, TH1* h_eff,
           const std::string& output_name)
    {
        TCanvas* c_combo = new TCanvas(("c_" + output_name).c_str(),
                                       "Events and Efficiency", 800, 600);

        TPad* pad1 = new TPad("pad1", "", 0, 0.37, 1, 1.0);
        pad1->SetTopMargin(0.08);
        pad1->SetBottomMargin(0.08);
        pad1->SetLeftMargin(0.15);
        pad1->SetRightMargin(0.08);
        pad1->Draw();

        TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.33);
        pad2->SetTopMargin(0.06);
        pad2->SetBottomMargin(0.15);
        pad2->SetLeftMargin(0.15);
        pad2->SetRightMargin(0.08);
        pad2->Draw();

        pad1->cd();
        h_true_pass->Draw("HIST E1");
        h_reco_pass->Draw("HIST E1 SAME");

        TLegend* leg = new TLegend(0.62, 0.75, 0.88, 0.88);
        leg->AddEntry(h_true_pass, "True passed events", "l");
        leg->AddEntry(h_reco_pass, "Reco passed events", "l");
        leg->Draw();

        pad2->cd();
        h_eff->Draw("HIST E1");

        c_combo->SaveAs(Form("Pi0Analysis/%s.png", output_name.c_str()));
    };

    // -------------------------------------------------------------------------
    // Draw matrices and combo plots (original behavior)
    // -------------------------------------------------------------------------
    DrawMatrix(h_response_E_idx, "pi0_response_matrix_energy");
    DrawMatrix(h_response_P_idx, "pi0_response_matrix_momentum");
    DrawMatrix(h_response_E,     "pi0_response_matrix_energy_unindexed");
    DrawMatrix(h_response_P,     "pi0_response_matrix_momentum_unindexed");

    // draw the MiniBooNE LEE-binned response matrix for inspection
    DrawMatrix(h_response_E_LEE, "pi0_response_matrix_energy_MiniBooNELEE");

    PlotEventsAndEfficiency(h_pass_true, h_pass_reco_energy,
                            h_eff_pi0,
                            "pi0_events_and_efficiency_energy");

    // Classification
    TCanvas *c_class = new TCanvas("c_class", "Pi0 Classification", 800, 800);
    h_pi0_nuance_class->SetFillColor(kOrange - 3);
    h_pi0_nuance_class->Draw("HIST");
    c_class->SaveAs("Pi0Analysis/Pi0_Source_Classification_MomT.png");

    // 1D distributions
    TCanvas* c_trueE = new TCanvas("c_trueE", "True Neutrino Energy", 800, 800);
    h_total_true->Draw("HIST");
    c_trueE->SaveAs("Pi0Analysis/Pi0_Total_True_Energy_MomT.png");

    TCanvas* c_trueMom = new TCanvas("c_trueMom", "True Pi0 Momentum", 800, 800);
    h_total_true_mom->Draw("HIST");
    c_trueMom->SaveAs("Pi0Analysis/Pi0_Total_True_Momentum_MomT.png");

    TCanvas* c_passE = new TCanvas("c_passE", "Passed Events by True Neutrino Energy", 800, 800);
    h_pass_true->Draw("HIST");
    c_passE->SaveAs("Pi0Analysis/Pi0_Passed_True_Energy_MomT.png");

    TCanvas* c_passMom = new TCanvas("c_passMom", "Passed Events by True Pi0 Momentum", 800, 800);
    h_pass_true_mom->Draw("HIST");
    c_passMom->SaveAs("Pi0Analysis/Pi0_Passed_True_Momentum_MomT.png");

    TCanvas* c_passRecoE = new TCanvas("c_passRecoE", "Passed Events by Reco Neutrino Energy", 800, 800);
    h_pass_reco_energy->Draw("HIST");
    c_passRecoE->SaveAs("Pi0Analysis/Pi0_Passed_Reco_Energy_MomT.png");

    // CC / NC response matrices
    TCanvas* c_ecc = new TCanvas("c_ecc", "Energy Response (CC)", 800, 800);
    h_response_E_CC->Draw("COLZ");
    c_ecc->SaveAs("Pi0Analysis/Pi0_Response_Energy_CC_MomT.png");

    TCanvas* c_enc = new TCanvas("c_enc", "Energy Response (NC)", 800, 800);
    h_response_E_NC->Draw("COLZ");
    c_enc->SaveAs("Pi0Analysis/Pi0_Response_Energy_NC_MomT.png");

    TCanvas* c_pcc = new TCanvas("c_pcc", "Momentum Response (CC)", 800, 800);
    h_response_P_CC->Draw("COLZ");
    h_response_P_CC->SetMarkerColor(kWhite);
    c_pcc->SaveAs("Pi0Analysis/Pi0_Response_Momentum_CC_MomT.png");

    TCanvas* c_pnc = new TCanvas("c_pnc", "Momentum Response (NC)", 800, 800);
    h_response_P_NC->Draw("COLZ");
    c_pnc->SaveAs("Pi0Analysis/Pi0_Response_Momentum_NC_MomT.png");

    // Ratio plot
    TCanvas* c_ratio = new TCanvas("c_ratio", "Pass RecoE / Total True Pi0 Momentum", 800, 800);
    h_ratio_passRecoE_totalTrueP->SetLineColor(kGreen + 2);
    h_ratio_passRecoE_totalTrueP->SetLineWidth(2);
    h_ratio_passRecoE_totalTrueP->Draw("HIST");
    c_ratio->SaveAs("Pi0Analysis/Pi0_Ratio_PassRecoE_to_TotalTrueP_MomT.png");

    // Efficiency vs true momentum / energy
    gStyle->SetOptStat(0);

    TCanvas* c_effMom = new TCanvas("c_effMom", "Efficiency vs True Pi0 Momentum", 800, 800);
    h_eff_pi0_mom->SetMinimum(0);
    c_effMom->SetLeftMargin(0.18);
    c_effMom->SetRightMargin(0.20);
    c_effMom->SetBottomMargin(0.18);
    c_effMom->SetTopMargin(0.1);
    h_eff_pi0_mom->SetLineColor(kBlue + 2);
    h_eff_pi0_mom->SetLineWidth(2);
    h_eff_pi0_mom->Draw("E1 P");
    c_effMom->SaveAs("Pi0Analysis/Pi0_Efficiency_True_Momentum_MomT.png");

    TCanvas* c_eff = new TCanvas("c_eff", "Efficiency vs True Pi0 Energy", 800, 800);
    h_eff_pi0->SetMinimum(0);
    c_eff->SetLeftMargin(0.18);
    c_eff->SetRightMargin(0.20);
    c_eff->SetBottomMargin(0.18);
    c_eff->SetTopMargin(0.1);
    h_eff_pi0->SetLineColor(kBlue + 2);
    h_eff_pi0->SetLineWidth(2);
    h_eff_pi0->Draw("E1 P");
    c_eff->SaveAs("Pi0Analysis/Pi0_Efficiency_True_Energy.png");

    TCanvas* c_effMom_CC = new TCanvas("c_effMom_CC", "Efficiency vs True Pi0 Momentum (CC)", 800, 800);
    c_effMom_CC->SetLeftMargin(0.18);
    c_effMom_CC->SetRightMargin(0.20);
    c_effMom_CC->SetBottomMargin(0.18);
    c_effMom_CC->SetTopMargin(0.1);
    h_eff_pi0_mom_CC->SetMinimum(0);
    h_eff_pi0_mom_CC->SetLineColor(kRed + 1);
    h_eff_pi0_mom_CC->SetLineWidth(2);
    h_eff_pi0_mom_CC->Draw("E1 P");
    c_effMom_CC->SaveAs("Pi0Analysis/Pi0_Efficiency_True_Momentum_CC_MomT.png");

    TCanvas* c_effMom_NC = new TCanvas("c_effMom_NC", "Efficiency vs True Pi0 Momentum (NC)", 800, 800);
    c_effMom_NC->SetLeftMargin(0.18);
    c_effMom_NC->SetRightMargin(0.20);
    c_effMom_NC->SetBottomMargin(0.18);
    c_effMom_NC->SetTopMargin(0.1);
    h_eff_pi0_mom_NC->SetMinimum(0);
    h_eff_pi0_mom_NC->SetLineColor(kBlue + 1);
    h_eff_pi0_mom_NC->SetLineWidth(2);
    h_eff_pi0_mom_NC->Draw("HIST");
    h_eff_pi0_mom_NC->Draw("E1 P");
    c_effMom_NC->SaveAs("Pi0Analysis/Pi0_Efficiency_True_Momentum_NC_MomT.png");

    TCanvas* c_event_counts = new TCanvas("c_event_counts", "Pi0 Event Category Counts", 800, 800);
    h_pi0_event_counts->SetFillColor(kAzure + 1);
    h_pi0_event_counts->SetLineColor(kBlack);
    h_pi0_event_counts->Draw("HIST TEXT");
    c_event_counts->SaveAs("Pi0Analysis/Pi0_Event_Counts.png");

    // Extra 1D distributions
    TCanvas* c_recoE_dist = new TCanvas("c_recoE_dist", "Reco Energy Distribution", 800, 800);
    h_pass_reco_energy->SetLineColor(kViolet + 1);
    h_pass_reco_energy->SetLineWidth(2);
    h_pass_reco_energy->Draw("HIST");
    c_recoE_dist->SaveAs("Pi0Analysis/Pi0_Reco_Energy_Distribution_MomT.png");

    TCanvas* c_passMom_dist = new TCanvas("c_passMom_dist", "Passed True Pi0 Momentum", 800, 800);
    h_pass_true_mom->SetLineColor(kGreen + 2);
    h_pass_true_mom->SetLineWidth(2);
    h_pass_true_mom->Draw("HIST");
    c_passMom_dist->SaveAs("Pi0Analysis/Pi0_Truth_Passed_Momentum_Distribution_MomT.png");

    TCanvas* c_totalMom_dist = new TCanvas("c_totalMom_dist", "Total True Pi0 Momentum", 800, 800);
    h_total_true_mom->SetLineColor(kBlue + 1);
    h_total_true_mom->SetLineWidth(2);
    h_total_true_mom->Draw("HIST");
    c_totalMom_dist->SaveAs("Pi0Analysis/Pi0_Truth_Total_Momentum_Distribution_MomT.png");

    // -------------------------------------------------------------------------
    // Print event counts
    // -------------------------------------------------------------------------
    std::cout << "------ Pi0 Event Category Counts ------" << std::endl;
    std::cout << "Total CC truth pi0s       : " << count_total_CC << std::endl;
    std::cout << "Total NC truth pi0s       : " << count_total_NC << std::endl;
    std::cout << "Passed CC pi0s (PassOsc)  : " << count_passed_CC << std::endl;
    std::cout << "Passed NC pi0s (PassOsc)  : " << count_passed_NC << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    // -------------------------------------------------------------------------
    // Write Technote-style 9-bin outputs to ROOT 
    // -------------------------------------------------------------------------
    TFile *fout = new TFile("Pi0_ResponseMatrices_MomT.root", "RECREATE");
    h_response_E->Write("Pi0_Response_Energy_MomT");
    h_response_P->Write("Pi0_Response_Momentum_MomT");
    h_total_true->Write("Pi0_Total_True_Energy_MomT");
    h_total_true_mom->Write("Pi0_Total_True_Momentum_MomT");
    h_pass_true->Write("Pi0_Passed_True_Energy_MomT");
    h_pass_true_mom->Write("Pi0_Passed_True_Momentum_MomT");
    h_eff_pi0_mom->Write("Pi0_Efficiency_True_Momentum_MomT");
    h_pass_reco_energy->Write("Pi0_Passed_Reco_Energy_MomT");
    h_ratio_passRecoE_totalTrueP->Write("Pi0_Ratio_PassRecoE_to_TotalTrueP_MomT");
    h_total_true_mom_CC->Write("Pi0_Total_True_Momentum_CC_MomT");
    h_total_true_mom_NC->Write("Pi0_Total_True_Momentum_NC_MomT");
    h_eff_pi0_mom_CC->Write("Pi0_Efficiency_True_Momentum_CC_MomT");
    h_eff_pi0_mom_NC->Write("Pi0_Efficiency_True_Momentum_NC_MomT");
    h_pi0_nuance_class->Write("Pi0_Classification_MomT");
    fout->Close();

    // -------------------------------------------------------------------------
    // Write MiniBooNE LEE–binned response matrix to SEPARATE ROOT file
    // -------------------------------------------------------------------------
    TFile *foutLEE = new TFile("Pi0_ResponseMatrix_MiniBooNELEE.root", "RECREATE");
    h_response_E_LEE->Write("Pi0_Response_Energy_MiniBooNELEE");
    foutLEE->Close();

    // -------------------------------------------------------------------------
    // Final summary
    // -------------------------------------------------------------------------
    std::cout << "Entries processed: " << totalEntries
              << ", pi0 events: " << pi0count << std::endl;
}
