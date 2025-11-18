#include "response_matrix_utils.h"

// ---------------------------------------
// Normalise a TH2D row-wise
// ---------------------------------------
void NormaliseByRow(TH2D* hist)
{
    for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
        double sum = 0.0;
        for (int ix = 1; ix <= hist->GetNbinsX(); ++ix)
            sum += hist->GetBinContent(ix, iy);

        if (sum > 0) {
            for (int ix = 1; ix <= hist->GetNbinsX(); ++ix)
                hist->SetBinContent(ix, iy, hist->GetBinContent(ix, iy) / sum);
        }
    }
}

// ---------------------------------------
// Copy contents and label with indices
// ---------------------------------------
void CopyAndLabelIndexed(TH2D* source, TH2D* target)
{
    int nx = source->GetNbinsX();
    int ny = source->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix)
        for (int iy = 1; iy <= ny; ++iy)
            target->SetBinContent(ix, iy, source->GetBinContent(ix, iy));

    for (int ix = 1; ix <= nx; ++ix)
        target->GetXaxis()->SetBinLabel(ix, std::to_string(ix).c_str());

    for (int iy = 1; iy <= ny; ++iy)
        target->GetYaxis()->SetBinLabel(iy, std::to_string(iy).c_str());
}

// ---------------------------------------
// Global ROOT style settings
// ---------------------------------------
void ConfigureGlobalStyle()
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kInvertedDarkBodyRadiator);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.15);
}

// ---------------------------------------
// Draw a 2D matrix
// ---------------------------------------
TCanvas* DrawMatrix2D(TH2* hist, const std::string& outDir,
                      const std::string& baseName, const char* ctitle,
                      bool showLabels)
{
    TCanvas* c = new TCanvas(baseName.c_str(), ctitle, 900, 800);
    if (showLabels) hist->LabelsOption("h");
    hist->Draw("COLZ");

    gSystem->mkdir(outDir.c_str(), kTRUE);
    std::string fname = outDir + "/" + baseName + ".png";
    c->SaveAs(fname.c_str());
    return c;
}

// ---------------------------------------
// Draw indexed 2D matrix
// ---------------------------------------
TCanvas* DrawIndexedMatrix2D(TH2* hist, const std::string& outDir,
                             const std::string& baseName, const char* ctitle)
{
    TCanvas* c = new TCanvas(baseName.c_str(), ctitle, 900, 800);
    hist->LabelsOption("h", "XY");
    hist->Draw("COLZ TEXT");

    gSystem->mkdir(outDir.c_str(), kTRUE);
    std::string fname = outDir + "/" + baseName + ".png";
    c->SaveAs(fname.c_str());
    return c;
}

// ---------------------------------------
// Draw 1D histogram
// ---------------------------------------
TCanvas* DrawHist1D(TH1* hist, const std::string& outDir,
                    const std::string& baseName, const char* ctitle,
                    Color_t color, int markerStyle,
                    const char* drawOpt,
                    double ymin, double ymax,
                    const char* title)
{
    TCanvas* c = new TCanvas(baseName.c_str(), ctitle, 800, 600);
    hist->SetLineColor(color);
    hist->SetMarkerStyle(markerStyle);

    if (ymin != ymax) { hist->SetMinimum(ymin); hist->SetMaximum(ymax); }
    if (title) hist->SetTitle(title);

    hist->Draw(drawOpt);

    gSystem->mkdir(outDir.c_str(), kTRUE);
    std::string fname = outDir + "/" + baseName + ".png";
    c->SaveAs(fname.c_str());
    return c;
}

// ---------------------------------------
// Combined event + efficiency panel
// ---------------------------------------
TCanvas* PlotEventsAndEfficiency(TH1* h_true_pass, TH1* h_reco_pass, TH1* h_eff,
                                 const std::string& outPath)
{
    TCanvas* c = new TCanvas("combo", "Events + Efficiency", 900, 900);

    TPad* pad1 = new TPad("pad1", "Events", 0, 0.35, 1, 1.0);
    pad1->SetBottomMargin(0.05);
    pad1->Draw();

    TPad* pad2 = new TPad("pad2", "Efficiency", 0, 0.0, 1, 0.32);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.20);
    pad2->Draw();

    // Top pad
    pad1->cd();
    h_true_pass->SetLineColor(kBlue + 1);
    h_true_pass->SetLineWidth(2);
    h_true_pass->Draw("HIST E1");

    h_reco_pass->SetLineColor(kRed + 1);
    h_reco_pass->SetLineWidth(2);
    h_reco_pass->Draw("HIST E1 SAME");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(h_true_pass, "True passed", "l");
    leg->AddEntry(h_reco_pass, "Reco passed", "l");
    leg->Draw();

    // Bottom pad
    pad2->cd();
    h_eff->SetMinimum(0.0);
    h_eff->SetMaximum(1.0);
    h_eff->SetLineColor(kBlack);
    h_eff->SetLineWidth(2);
    h_eff->Draw("HIST E1");

    c->SaveAs(outPath.c_str());
    return c;
}
