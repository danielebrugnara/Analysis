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
#include <TCanvas.h>


class SpectrumAnalyzer {
public:
    explicit SpectrumAnalyzer(const std::string &, const bool&, const std::string&, const bool&);

    ~SpectrumAnalyzer();

    void Analyze();

private:
    typedef std::vector<std::pair<double, double>> IntensityData;

    std::string file_name;
    bool debug_canvas;
    bool use_main_transitions;
    bool simulation;
    TH2D gg;
    TH1D hspec;
    std::vector<TH1D> cry_spec;
    TVector2 start_stop;
    double nevts;
    Plotter plotter;

    TGraphErrors relative_eff_graph;
    std::vector<TGraphErrors*> cry_relative_eff_graph;
    TGraphErrors absolute_eff_graph;
    TGraphErrors sigma_graph;
    TGraphErrors tau_graph;

    LevelScheme levelscheme;
    IntensityData eu152_intensities;
    std::vector<std::pair<const Gamma *, const Gamma *>> gamma_gamma;
    const double fit_interval;
    const double proj_interval;
    unsigned fit_counter;

    double GetPeakIntegral(TH1D &, const double &);
    static IntensityData ReadIntensities(const std::string&);
    void GenerateRelativeEffGraph(const TH1D&, TGraphErrors&, TGraphErrors&, TGraphErrors&, int);
    void GenerateAbsoluteEffGraph();
    void GenerateCrystalEffGraph();

    static TF1 FitEffCurve(TGraphErrors&);

    static void UpdateErrors(TH1D &);

    static TH1D* ProjectAndSubtract(const TH2D &,
                             const double &, const double &,
                             const double &, const double &,
                             const double &, const double &);

};
#endif //EFFANALYSIS_SPECTRUMANALYZER_H
