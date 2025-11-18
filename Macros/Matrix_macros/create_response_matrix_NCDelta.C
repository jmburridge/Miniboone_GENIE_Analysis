#include "response_matrix_utils.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

void create_response_matrix_NCDelta()
{
    ConfigureGlobalStyle();

    // -----------------------
    // Input file
    // -----------------------
    TFile* file = TFile::Open("../MiniBooNEDatasets2023/miniboone_mc_all_ncdelta.root");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open NC delta input file!" << std::endl;
        return;
    }

    TTree* tree = (TTree*)file->Get("MiniBooNE_CCQE");
    if (!tree) {
        std::cerr << "ERROR: TTree MiniBooNE_CCQE missing!" << std::endl;
        return;
    }

    // -----------------------
    // Native NCΔ binning
    // -----------------------
    const int nbins_reco = 8;
    double reco_bins[9] = {
        0.200, 0.300, 0.375, 0.475, 0.550,
        0.675, 0.800, 1.000, 1.200
    };

    const int nbins_true = 7;
    double true_bins[8] = {
        0.250, 0.500, 0.775, 1.000,
        1.275, 1.500, 2.000, 3.000
    };

    // LEE binning
    const int nbinsLEE = 11;
    double lee_bins[12] = {
        0.2, 0.3, 0.375, 0.475, 0.55,
        0.675, 0.8, 0.95, 1.1, 1.3,
        1.5, 3.0
    };

    // -----------------------
    // Histograms
    // -----------------------
    TH1D* h_true = new TH1D("h_true", "True E", nbins_true, true_bins);
    TH1D* h_pass = new TH1D("h_pass", "Passed True E", nbins_true, true_bins);
    TH1D* h_reco = new TH1D("h_reco", "Reco E", nbins_reco, reco_bins);

    TH2D* h_smear = new TH2D("h_smear", "Smearing", nbins_reco, reco_bins,
                                                      nbins_true, true_bins);
    TH2D* h_resp  = new TH2D("h_resp",  "Response", nbins_reco, reco_bins,
                                                      nbins_true, true_bins);

    // LEE
    TH1D* h_true_LEE = new TH1D("h_true_LEE", "True E LEE", nbinsLEE, lee_bins);
    TH1D* h_pass_LEE = new TH1D("h_pass_LEE", "Passed True E LEE", nbinsLEE, lee_bins);
    TH1D* h_reco_LEE = new TH1D("h_reco_LEE", "Reco E LEE", nbinsLEE, lee_bins);
    TH2D* h_smear_LEE = new TH2D("h_smear_LEE", "Smearing LEE",
            nbinsLEE, lee_bins, nbinsLEE, lee_bins);
    TH2D* h_resp_LEE  = new TH2D("h_resp_LEE", "Response LEE",
            nbinsLEE, lee_bins, nbinsLEE, lee_bins);

    // -----------------------
    // Tree variables
    // -----------------------
    float NuMomT = 0;
    float RecoEnuQE = 0;
    bool PassOsc = false;
    float Weight = 1;

    tree->SetBranchAddress("NuMomT", &NuMomT);
    tree->SetBranchAddress("RecoEnuQE", &RecoEnuQE);
    tree->SetBranchAddress("PassOsc", &PassOsc);
    if (tree->GetBranch("Weight"))
        tree->SetBranchAddress("Weight", &Weight);

    // -----------------------
    // Event loop
    // -----------------------
    Long64_t N = tree->GetEntries();
    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        h_true->Fill(NuMomT, Weight);
        h_true_LEE->Fill(NuMomT, Weight);

        if (PassOsc) {
            h_pass->Fill(NuMomT, Weight);
            h_reco->Fill(RecoEnuQE, Weight);
            h_smear->Fill(RecoEnuQE, NuMomT, Weight);

            h_pass_LEE->Fill(NuMomT, Weight);
            h_reco_LEE->Fill(RecoEnuQE, Weight);
            h_smear_LEE->Fill(RecoEnuQE, NuMomT, Weight);
        }
    }

    // -----------------------
    // Normalisation
    // -----------------------
    NormaliseByRow(h_smear);
    NormaliseByRow(h_smear_LEE);

    h_resp->Add(h_smear);
    h_resp_LEE->Add(h_smear_LEE);

    // -----------------------
    // Output files
    // -----------------------
    TFile* fout = new TFile("response_ncdeltas.root", "RECREATE");
    h_true->Write();
    h_pass->Write();
    h_reco->Write();
    h_smear->Write();
    h_resp->Write();
    fout->Close();

    TFile* foutLEE = new TFile("response_ncdeltas_LEE.root", "RECREATE");
    h_true_LEE->Write();
    h_pass_LEE->Write();
    h_reco_LEE->Write();
    h_smear_LEE->Write();
    h_resp_LEE->Write();
    foutLEE->Close();

    std::cout << "Finished NC DELTAS response matrices" << std::endl;
}
