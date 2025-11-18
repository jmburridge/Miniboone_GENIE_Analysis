#ifndef RESPONSE_MATRIX_UTILS_H
#define RESPONSE_MATRIX_UTILS_H

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TSystem.h"
#include <string>

// ------------------------------
// Normalisation & utility helpers
// ------------------------------
void NormaliseByRow(TH2D* hist);
void CopyAndLabelIndexed(TH2D* source, TH2D* target);

// ------------------------------
// Style configuration
// ------------------------------
void ConfigureGlobalStyle();

// ------------------------------
// Drawing helpers
// ------------------------------
TCanvas* DrawMatrix2D(TH2* hist, const std::string& outDir,
                      const std::string& baseName, const char* ctitle,
                      bool showLabels = false);

TCanvas* DrawIndexedMatrix2D(TH2* hist, const std::string& outDir,
                             const std::string& baseName, const char* ctitle);

TCanvas* DrawHist1D(TH1* hist, const std::string& outDir,
                    const std::string& baseName, const char* ctitle,
                    Color_t color, int markerStyle,
                    const char* drawOpt = "HIST E1",
                    double ymin = -1, double ymax = -1,
                    const char* title = nullptr);

TCanvas* PlotEventsAndEfficiency(TH1* h_true_pass, TH1* h_reco_pass, TH1* h_eff,
                                 const std::string& outPath);

#endif
