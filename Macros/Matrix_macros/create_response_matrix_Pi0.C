// ============================================================================
// Final Clean π0 Response Matrix Macro (Improved)
// Energy-only, no CC/NC, no momentum, no NUANCE classification
// Standardised ROOT outputs only:
//   - response_pi0.root (9-bin technote format)
//   - response_pi0_LEE.root (11-bin LEE format)
// ROOT files saved to:
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Pi0/Pi0_Root_Files/
// PNGs saved to:
//   /exp/uboone/app/users/jburridg/Geometry/Analysis/Pi0/Pi0_Histograms/
// Professional PNG naming + original colour scheme + stat boxes + improved margins
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

// -----------------------------------------------------------------------------
void create_response_matrix_Pi0()
{
    // -----------------------------------------------------------------------------
    // GLOBAL STYLE CONFIGURATION (colour, stats, margins)
    // -----------------------------------------------------------------------------
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    gStyle->SetOptStat(1110);      // entries, mean, RMS
    gStyle->SetStatBorderSize(1);
    gStyle->SetStatX(0.88);
    gStyle->SetStatY(0.88);

    // -----------------------------------------------------------------------------
    // Binning
    // -----------------------------------------------------------------------------
    const int nbins_true = 10;
    const int nbins_reco = 9;
    double recoE_bins[nbins_rec0+1] = {0.150,0.200,0.300,0.375,0.475,0.550,0.675,0.800,1.000,1.200};
    double trueE_bins[nbins_true+1] = {0.250,0.500,0.600,0.700,0.900,1.100,1.250,1.500,2.0,3.0};

    const int lee_nbins = 11;
    double lee_bins[lee_nbins+1] = {0.2,0.3,0.375,0.475,0.55,0.675,0.8,0.95,1.1,1.3,1.5,3.0};

    // -----------------------------------------------------------------------------
    // Energy-only histograms
    // -----------------------------------------------------------------------------
    TH1D* h_total_true = new TH1D("h_total_true",
        "Total True E_{#nu};E_{#nu} [GeV];Events", nbins, trueE_bins);

    TH1D* h_pass_true = new TH1D("h_pass_true",
        "Passed True E_{#nu};E_{#nu} [GeV];Events", nbins, trueE_bins);

    TH1D* h_pass_reco = new TH1D("h_pass_reco",
        "Passed Reco E_{#nu}^{QE};Reco E_{#nu}^{QE} [GeV];Events",
        nbins, recoE_bins);

    TH2D* h_response_E = new TH2D("h_response_E",
        "Response Matrix (Reco vs True E);Reco E_{#nu}^{QE};True E_{#nu}",
        nbins, recoE_bins, nbins, trueE_bins);

    TH1D* h_eff_pi0 = new TH1D("h_eff_pi0",
        "Pi0 Efficiency;True E_{#nu} [GeV];Efficiency",
        nbins, trueE_bins);

    // -----------------------------------------------------------------------------
    // LEE-binned histograms
    // -----------------------------------------------------------------------------
    TH1D* h_total_true_LEE = new TH1D("h_total_true_LEE",
        "Total True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);

    TH1D* h_pass_true_LEE = new TH1D("h_pass_true_LEE",
        "Passed True E_{#nu} (LEE);E_{#nu} [GeV];Events", lee_nbins, lee_bins);

    TH2D* h_response_E_LEE = new TH2D("h_response_E_LEE",
        "LEE Response Matrix;Reco E_{#nu}^{QE};True E_{#nu}",
        lee_nbins, lee_bins, lee_nbins, lee_bins);

    long long totalEntries = 0, pi0count = 0;

    // -----------------------------------------------------------------------------
    // File loop
    // -----------------------------------------------------------------------------
    for (int fileIndex=1; fileIndex<=10; fileIndex++)
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

        std::vector<int> *FSPType=nullptr;
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

        for (int i=0; i<N; i++)
        {
            t->GetEntry(i);
            if (!CheckPointers(FSPType,Vx,Vy,Vz,MomX,MomY,MomZ,MomT))
                continue;

            unsigned npi0 = sp::Pi0Details(
                NFSP,*FSPType,*Vx,*Vy,*Vz,*MomX,*MomY,*MomZ,*MomT);

            if (npi0==0) continue;

            StackedBkgdType_t bkgd =
                StackHistoBkgd(false,true,(NuanceType_t)NUANCEChan,
                               (NuType_t)NuType,(GEANT3Type_t)NuParentID);

            if (bkgd != kBKGD_PI0) continue;

            pi0count++;

            h_total_true->Fill(NuMomT,Weight);
            h_total_true_LEE->Fill(NuMomT,Weight);

            if (PassOsc) {
                h_pass_true->Fill(NuMomT,Weight);
                h_pass_true_LEE->Fill(NuMomT,Weight);

                h_pass_reco->Fill(RecoEnuQE,Weight);

                h_response_E->Fill(RecoEnuQE,NuMomT,Weight);
                h_response_E_LEE->Fill(RecoEnuQE,NuMomT,Weight);
            }
        }

        f->Close();
    }

    // -----------------------------------------------------------------------------
    // Normalise response matrices
    // -----------------------------------------------------------------------------
    for (int iy=1; iy<nbins+1; iy++) {
        double den = h_total_true->GetBinContent(iy);
        if (den<=0) continue;
        for (int ix=1; ix<nbins+1; ix++)
            h_response_E->SetBinContent(ix,iy,
                h_response_E->GetBinContent(ix,iy)/den);
    }

    for (int iy=1; iy<=lee_nbins; iy++) {
        double den = h_total_true_LEE->GetBinContent(iy);
        if (den<=0) continue;
        for (int ix=1; ix<=lee_nbins; ix++)
            h_response_E_LEE->SetBinContent(ix,iy,
                h_response_E_LEE->GetBinContent(ix,iy)/den);
    }

    // -----------------------------------------------------------------------------
    // Efficiency
    // -----------------------------------------------------------------------------
    for (int i=1; i<=nbins; i++) {
        double num = h_pass_true->GetBinContent(i);
        double den = h_total_true->GetBinContent(i);
        if (den>0) {
            double eff = num/den;
            double err = sqrt(eff*(1-eff)/den);
            h_eff_pi0->SetBinContent(i,eff);
            h_eff_pi0->SetBinError(i,err);
        }
    }

    // -----------------------------------------------------------------------------
    // Indexed response matrix
    // -----------------------------------------------------------------------------
    TH2D* h_response_E_idx =
        new TH2D("h_response_E_idx",";Reco E bin;True E bin",
                 nbins,0,nbins, nbins,0,nbins);
    CopyAndLabelIndexed(h_response_E,h_response_E_idx);

    // -----------------------------------------------------------------------------
    // Output directories
    // -----------------------------------------------------------------------------
    const char* png_dir =
        "/exp/uboone/app/users/jburridg/Geometry/Analysis/Pi0/Pi0_Histograms/";
    const char* root_dir =
        "/exp/uboone/app/users/jburridg/Geometry/Analysis/Pi0/Pi0_Root_Files/";

    // -----------------------------------------------------------------------------
    // Plotting utilities (improved margins)
    // -----------------------------------------------------------------------------
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
        TCanvas* c=new TCanvas(fname,fname,800,800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.15);
        c->SetBottomMargin(0.14);
        c->SetTopMargin(0.08);
        h->Draw("HIST E");
        std::string out = std::string(png_dir) + fname + ".png";
        c->SaveAs(out.c_str());
    };

    // -----------------------------------------------------------------------------
    // Save PNGs
    // -----------------------------------------------------------------------------
    SaveMatrix(h_response_E,     "Pi0_ResponseMatrix_Energy");
    SaveMatrix(h_response_E_idx, "Pi0_ResponseMatrix_Energy_Indexed");
    SaveMatrix(h_response_E_LEE,"Pi0_ResponseMatrix_Energy_LEE");

    Save1D(h_total_true, "Pi0_TrueEnergy_Total");
    Save1D(h_pass_true,  "Pi0_TrueEnergy_Passed");
    Save1D(h_pass_reco,  "Pi0_RecoEnergy_Passed");
    Save1D(h_eff_pi0,    "Pi0_Efficiency_TrueEnergy");

    // -----------------------------------------------------------------------------
    // Standardised export (9-bin)
    // -----------------------------------------------------------------------------
    {
        TH1D* h_true_export = (TH1D*)h_total_true->Clone("h_true");
        TH1D* h_pass_export = (TH1D*)h_pass_true->Clone("h_pass");
        TH1D* h_reco_export = (TH1D*)h_pass_reco->Clone("h_reco");

        TH2D* h_smear_unorm = (TH2D*)h_response_E->Clone("h_smear_unorm");
        TH2D* h_resp        = (TH2D*)h_response_E->Clone("h_resp");

        TH2D* h_resp_index =
            new TH2D("h_resp_index","Response indexed",
                     nbins,0,nbins, nbins,0,nbins);
        CopyAndLabelIndexed(h_resp,h_resp_index);

        std::string outfile = std::string(root_dir) + "response_pi0.root";
        TFile out(outfile.c_str(),"RECREATE");
        h_true_export->Write();
        h_pass_export->Write();
        h_reco_export->Write();
        h_smear_unorm->Write();
        h_resp->Write();
        h_resp_index->Write();
        out.Close();
    }

    // -----------------------------------------------------------------------------
    // Standardised export (LEE, 11-bin)
    // -----------------------------------------------------------------------------
    {
        TH1D* h_true_LEE_export = (TH1D*)h_total_true_LEE->Clone("h_true_LEE");
        TH1D* h_pass_LEE_export = (TH1D*)h_pass_true_LEE->Clone("h_pass_LEE");

        TH1D* h_reco_LEE = new TH1D("h_reco_LEE","Reco LEE",lee_nbins,lee_bins);
        for (int ix=1; ix<=lee_nbins; ix++){
            double sum=0;
            for (int iy=1; iy<=lee_nbins; iy++)
                sum += h_response_E_LEE->GetBinContent(ix,iy);
            h_reco_LEE->SetBinContent(ix,sum);
        }

        TH2D* h_smear_LEE_unorm = (TH2D*)h_response_E_LEE->Clone("h_smear_LEE_unorm");
        TH2D* h_resp_LEE        = (TH2D*)h_response_E_LEE->Clone("h_resp_LEE");

        TH2D* h_resp_index_LEE =
            new TH2D("h_resp_index_LEE","LEE Response indexed",
                     lee_nbins,0,lee_nbins, lee_nbins,0,lee_nbins);
        CopyAndLabelIndexed(h_resp_LEE,h_resp_index_LEE);

        std::string outfile = std::string(root_dir) + "response_pi0_LEE.root";
        TFile out(outfile.c_str(),"RECREATE");
        h_true_LEE_export->Write();
        h_pass_LEE_export->Write();
        h_reco_LEE->Write();
        h_smear_LEE_unorm->Write();
        h_resp_LEE->Write();
        h_resp_index_LEE->Write();
        out.Close();
    }

    // -----------------------------------------------------------------------------
    std::cout
        << "[PI0] Final energy-only π0 response matrices completed.\n"
        << "      Entries processed: " << totalEntries << "\n"
        << "      π0 events kept:    " << pi0count << std::endl;
}
