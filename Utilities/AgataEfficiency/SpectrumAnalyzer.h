//
// Created by daniele on 13/07/2020.
//

#ifndef EFFANALYSIS_SPECTRUMANALYZER_H
#define EFFANALYSIS_SPECTRUMANALYZER_H

#include <string>

#include "DiaGraph.h"
#include "LevelScheme.h"
#include "Plotter.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFitResult.h>


class SpectrumAnalyzer {
public:
    explicit SpectrumAnalyzer(const std::string &);

    ~SpectrumAnalyzer();

    void Analyze();

private:
    typedef std::vector<std::pair<double, double>> IntensityData;

    TH2D gg;
    Plotter plotter;
    TGraph effgraph;
    TGraph relative_effgraph;
    //void SetupLevelSchemes();
    LevelScheme levelscheme;
    IntensityData eu152_intensities;
    std::vector<std::pair<const Gamma *, const Gamma *>> gamma_gamma;
    const double fit_interval;
    const double proj_interval;
    unsigned fit_counter;

    double GetPeakIntegral(TH1D &, const double &);
    static IntensityData ReadIntensities(const std::string&);
    TGraph GenerateRelativeEffGraph();

    static void UpdateErrors(TH1D &);

    static TH1D* ProjectAndSubtract(const TH2D &,
                             const double &, const double &,
                             const double &, const double &,
                             const double &, const double &);
    static TH1D* Subtract(TH1D&);

};
#endif //EFFANALYSIS_SPECTRUMANALYZER_H
