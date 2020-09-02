//
// Created by daniele on 13/07/2020.
//

#ifndef EFFANALYSIS_SPECTRUMANALYZER_H
#define EFFANALYSIS_SPECTRUMANALYZER_H

#include <string>
#include <algorithm>
#include <cmath>

#include "DiaGraph.h"
#include "LevelScheme.h"
#include "Plotter.h"
#include "Fitter.h"

#include <TFile.h>
#include <TVirtualFitter.h>
#include <TCanvas.h>
#include <Math/MinimizerOptions.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TVector2.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TFitResult.h>


class SpectrumAnalyzer {
public:
    explicit SpectrumAnalyzer(const std::string &, const bool&);

    ~SpectrumAnalyzer();

    void Analyze();

private:
    typedef std::vector<std::pair<double, double>> IntensityData;

    bool debug_canvas;
    bool simulation;
    TH2D gg;
    TH1D hspec;
    TVector2 start_stop;
    double nevts;
    Plotter plotter;
    TGraph effgraph;
    TGraphErrors sigmagraph;
    TGraphErrors sigmagaussgraph;
    TGraphErrors relative_effgraph;
    TGraphErrors scaled_relative_effgraph;
    TGraphErrors relative_integralgraph;
    TGraph relative_intgraph;
    //void SetupLevelSchemes();
    LevelScheme levelscheme;
    IntensityData eu152_intensities;
    std::vector<std::pair<const Gamma *, const Gamma *>> gamma_gamma;
    const double fit_interval;
    const double proj_interval;
    unsigned fit_counter;

    double GetPeakIntegral(TH1D &, const double &);
    static IntensityData ReadIntensities(const std::string&);
    void GenerateRelativeEffGraph();

    static void UpdateErrors(TH1D &);

    static TH1D* ProjectAndSubtract(const TH2D &,
                             const double &, const double &,
                             const double &, const double &,
                             const double &, const double &);
    static TH1D* Subtract(TH1D&);

};
#endif //EFFANALYSIS_SPECTRUMANALYZER_H