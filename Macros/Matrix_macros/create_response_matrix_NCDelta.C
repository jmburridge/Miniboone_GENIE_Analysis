// ============================================================================
// Final Clean NCΔ Response Matrix Macro (Parallel to π0 & Nue Macros)
// Energy-only response matrices
// Standardised ROOT outputs only:
//     • response_ncdelta.root        (NCΔ technote binning)
//     • response_ncdelta_LEE.root    (MiniBooNE LEE binning)
// ROOT files saved to:
//     /exp/uboone/app/users/jburridg/Geometry/Analysis/NCDelta/NCDelta_Root_Files/
// PNGs saved to:
//     /exp/uboone/app/users/jburridg/Geometry/Analysis/NCDelta/NCDelta_Histograms/
// ----------------------------------------------------------------------------
// Signal definition:
//     StackedBkgdType_t == kBKGD_DELTA
// ----------------------------------------------------------------------------
// Style: original colour scheme, stat boxes ON, improved margins,
// professional PNG naming.
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

#include "../../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedTypes.h"
#include "../../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.h"
#include "../../MiniBooNEDatasets2023/CombinedFunctions_from_Fortran/CombinedFunctions.cxx"

using namespace sp;

// ----------------------------------------------------------------------------
// Utility functions
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// Main macro
// ----------------------------------------------------------------------------
void create_response_matrix_NCDelta()
{
    // ----------------------------------------------------------------------------
    // Global style
    // ----------------------------------------------------------------------------
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    gStyle->SetOptStat(1110);
    gStyle->SetStatX(0.88);
    gStyle->SetStatY(0.88);
    gStyle->SetStatBorderSize(1);

    // ----------------------------------------------------------------------------
    // NCΔ Binning (technote)
    // ----------------------------------------------------------------------------
    const int nbins_reco = 8;
    double reco_bins[9] =
        {0.200, 0.300, 0.375, 0.475, 0.550, 0.675, 0.800, 1.000, 1.200};

    const int nbins_true = 7;
    double true_bins[8] =
        {0.250, 0.500, 0.775, 1.000, 1.275, 1.500, 2.000, 3.000};

    // ----------------------------------------------------------------------------
    // MiniBooNE LEE binning (same as pi0/nue)
    // ----------------------------------------------------------------------------
    const int nbinsLEE = 11;
    double lee_bins[12] =
        {0.2,0.3,0.375,0.475,0.55,0.675,0.8,0.95,1.1,1.3,1.5,3.0};

    // ----------------------------------------------------------------------------
    // Histograms (technote)
    // ----------------------------------------------------------------------------
    TH1D* h_true   = new TH1D("h_true",
                    "Total True E_{#nu};E_{#nu} [GeV];Events",
                    nbins_true, true_bins);

    TH1D* h_pass   = new TH1D("h_pass",
                    "Passed True E_{#nu};E_{#nu} [GeV];Events",
                    nbins_true, true_bins);

    TH1D* h_reco   = new TH1D("h_reco",
                    "Passed Reco E_{#nu}^{QE};Reco E_{#nu}^{QE} [GeV];Events",
                    nbins_reco, reco_bins);

    TH2D* h_smear  = new TH2D("h_smear",
                    "Raw Smearing Matrix;Reco E_{#nu}^{QE};True E_{#nu}",
                    nbins_reco, reco_bins,
                    nbins_true, true_bins);

    TH2D* h_resp   = new TH2D("h_resp",
                    "Response Matrix (Reco vs True E);Reco E_{#nu}^{QE};True E_{#nu}",
                    nbins_reco, reco_bins,
                    nbins_true, true_bins);

    // ----------------------------------------------------------------------------
    // LEE histograms
    // ----------------------------------------------------------------------------
    TH1D* h_true_LEE  = new TH1D("h_true_LEE",
                       "Total True E_{#nu} (LEE);E_{#nu} [GeV];Events",
                       nbinsLEE, lee_bins);

    TH1D* h_pass_LEE  = new TH1D("h_pass_LEE",
                       "Passed True E_{#nu} (LEE);E_{#nu} [GeV];Events",
                       nbinsLEE, lee_bins);

    TH2D* h_smear_LEE = new TH2D("h_smear_LEE",
                       "Raw Smearing (LEE);Reco E_{#nu}^{QE};True E_{#nu}",
                       nbinsLEE, lee_bins,
                       nbinsLEE, lee_bins);

    TH2D* h_resp_LEE  = new TH2D("h_resp_LEE",
                       "Response (LEE);Reco E_{#nu}^{QE};True E_{#nu}",
                       nbinsLEE, lee_bins,
                       nbinsLEE, lee_bins);

    // ----------------------------------------------------------------------------
    long long totalEntries = 0, kept = 0;

    // ----------------------------------------------------------------------------
    // File loop
    // ----------------------------------------------------------------------------
    for (int fileIndex = 1; fileIndex <= 10; fileIndex++)
    {
        std::stringstream ss;
        ss << "../../MiniBooNEDatasets2023/output_osc_mc_detail_"
           << fileIndex << ".root";

        TFile* f = TFile::Open(ss.str().c_str());
        if (!f || f->IsZombie()) continue;

        TTree* t = (TTree*)f->Get("MiniBooNE_CCQE");
        if (!t) { f->Close(); continue; }

        int NFSP, NUANCEChan, NuType, NuParentID;
        float Energy, RecoEnuQE, Weight, NuMomT;
        bool PassOsc;

        std::vector<int>* FSPType = nullptr;
        std::vector<float> *Vx=nullptr, *Vy=nullptr, *Vz=nullptr;
        std::vector<float> *MomX=nullptr, *MomY=nullptr, *MomZ=nullptr, *MomT=nullptr;

        // Branches
        t->SetBranchAddress("NFSP",&NFSP);
        t->SetBranchAddress("FSPType",&FSPType);
        t->SetBranchAddress("VertexX",&Vx);
        t->SetBranchAddress("VertexY",&Vy);
        t->SetBranchAddress("VertexZ",&Vz);
        t->SetBranchAddress("MomX",&MomX);
        t->SetBranchAddress("MomY",&MomY);
        t->SetBranchAddress("MomZ",&MomZ);
        t->SetBranchAddress("MomT",&MomT);
        t->SetBranchAddress("NUANCEChan",&NUANCEChan);
        t->SetBranchAddress("NuType",&NuType);
        t->SetBranchAddress("NuMomT",&NuMomT);
        t->SetBranchAddress("NuParentID",&NuParentID);
        t->SetBranchAddress("Energy",&Energy);
        t->SetBranchAddress("RecoEnuQE",&RecoEnuQE);
        t->SetBranchAddress("Weight",&Weight);
        t->SetBranchAddress("PassOsc",&PassOsc);

        int N = t->GetEntries();
        totalEntries += N;

        for (int i = 0; i < N; i++)
        {
            t->GetEntry(i);
            if (!CheckPointers(FSPType,Vx,Vy,Vz,MomX,MomY,MomZ,MomT))
                continue;

            // Select NCΔ
            StackedBkgdType_t bkgd =
                StackHistoBkgd(false,false,(NuanceType_t)NUANCEChan,
                               (NuType_t)NuType,(GEANT3Type_t)NuParentID);

            if (bkgd != kBKGD_DELTA) continue;

            kept++;

            // Fill truth
            h_true->Fill(NuMomT,Weight);
            h_true_LEE->Fill(NuMomT,Weight);

            if (PassOsc)
            {
                h_pass->Fill(NuMomT,Weight);
                h_pass_LEE->Fill(NuMomT,Weight);

                h_reco->Fill(RecoEnuQE,Weight);

                h_smear->Fill(RecoEnuQE,NuMomT,Weight);
                h_smear_LEE->Fill(RecoEnuQE,NuMomT,Weight);
            }
        }

        f->Close();
    }

    // ----------------------------------------------------------------------------
    // Column-normalise response matrices
    // ----------------------------------------------------------------------------
    for (int iy=1; iy<=nbins_true; iy++)
    {
        double den = h_true->GetBinContent(iy);
        if (den <= 0) continue;

        for (int ix=1; ix<=nbins_reco; ix++)
            h_resp->SetBinContent(
                ix, iy,
                h_smear->GetBinContent(ix,iy) / den
            );
    }

    for (int iy=1; iy<=nbinsLEE; iy++)
    {
        double den = h_true_LEE->GetBinContent(iy);
        if (den <= 0) continue;

        for (int ix=1; ix<=nbinsLEE; ix++)
        {
            h_resp_LEE->SetBinContent(
                ix, iy,
                h_smear_LEE->GetBinContent(ix,iy) / den
            );
        }
    }

    // ----------------------------------------------------------------------------
    // Indexed matrices
    // ----------------------------------------------------------------------------
    TH2D* h_resp_index =
        new TH2D("h_resp_index",
                 ";Reco bin;True bin",
                 nbins_reco,0,nbins_reco,
                 nbins_true,0,nbins_true);
    CopyAndLabelIndexed(h_resp,h_resp_index);

    TH2D* h_resp_index_LEE =
        new TH2D("h_resp_index_LEE",
                 ";Reco bin;True bin",
                 nbinsLEE,0,nbinsLEE,
                 nbinsLEE,0,nbinsLEE);
    CopyAndLabelIndexed(h_resp_LEE,h_resp_index_LEE);

    // ----------------------------------------------------------------------------
    // Output directories
    // ----------------------------------------------------------------------------
    const char* png_dir =
        "/exp/uboone/app/users/jburridg/Geometry/Analysis/NCDelta/NCDelta_Histograms/";
    const char* root_dir =
        "/exp/uboone/app/users/jburridg/Geometry/Analysis/NCDelta/NCDelta_Root_Files/";

    // ----------------------------------------------------------------------------
    // Plotting utilities
    // ----------------------------------------------------------------------------
    auto SaveMatrix = [&](TH2* h, const char* fname){
        TCanvas* c = new TCanvas(fname,fname,900,800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.18);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.08);
        h->LabelsOption("h");
        h->Draw("COLZ");
        std::string out = std::string(png_dir) + fname + ".png";
        c->SaveAs(out.c_str());
    };

    auto Save1D = [&](TH1* h, const char* fname){
        TCanvas* c = new TCanvas(fname,fname,800,800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.15);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.08);
        h->Draw("HIST E");
        std::string out = std::string(png_dir) + fname + ".png";
        c->SaveAs(out.c_str());
    };

    // ----------------------------------------------------------------------------
    // Save PNGs
    // ----------------------------------------------------------------------------
    SaveMatrix(h_resp,        "NCDelta_ResponseMatrix_Energy");
    SaveMatrix(h_resp_index,  "NCDelta_ResponseMatrix_Energy_Indexed");
    SaveMatrix(h_resp_LEE,    "NCDelta_ResponseMatrix_Energy_LEE");
    SaveMatrix(h_resp_index_LEE,"NCDelta_ResponseMatrix_Energy_LEE_Indexed");

    Save1D(h_true,     "NCDelta_TrueEnergy_Total");
    Save1D(h_pass,     "NCDelta_TrueEnergy_Passed");
    Save1D(h_reco,     "NCDelta_RecoEnergy_Passed");

    // ----------------------------------------------------------------------------
    // Standardised export: NCΔ
    // ----------------------------------------------------------------------------
    {
        std::string outfile = std::string(root_dir) + "response_ncdelta.root";
        TFile fexp(outfile.c_str(),"RECREATE");

        h_true->Write("h_true");
        h_pass->Write("h_pass");
        h_reco->Write("h_reco");
        h_smear->Write("h_smear_unorm");
        h_resp->Write("h_resp");
        h_resp_index->Write("h_resp_index");

        fexp.Close();
    }

    // ----------------------------------------------------------------------------
    // Standardised export: NCΔ (LEE)
    // ----------------------------------------------------------------------------
    {
        // Build reco LEE
        TH1D* h_reco_LEE = new TH1D("h_reco_LEE","Reco LEE",nbinsLEE,lee_bins);
        for (int ix=1; ix<=nbinsLEE; ix++)
        {
            double sum = 0;
            for (int iy=1; iy<=nbinsLEE; iy++)
                sum += h_smear_LEE->GetBinContent(ix,iy);
            h_reco_LEE->SetBinContent(ix,sum);
        }

        std::string outfile = std::string(root_dir) + "response_ncdelta_LEE.root";
        TFile fexp(outfile.c_str(),"RECREATE");

        h_true_LEE->Write("h_true_LEE");
        h_pass_LEE->Write("h_pass_LEE");
        h_reco_LEE->Write("h_reco_LEE");
        h_smear_LEE->Write("h_smear_LEE_unorm");
        h_resp_LEE->Write("h_resp_LEE");
        h_resp_index_LEE->Write("h_resp_index_LEE");

        fexp.Close();
    }

    // ----------------------------------------------------------------------------
    std::cout << "[NCDelta] Response matrices complete.\n"
              << "           Entries processed: " << totalEntries << "\n"
              << "           NCΔ events kept:   " << kept << std::endl;
}
