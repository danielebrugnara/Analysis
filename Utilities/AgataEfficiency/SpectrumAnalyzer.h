//
// Created by daniele on 13/07/2020.
//

#ifndef EFFANALYSIS_SPECTRUMANALYZER_H
#define EFFANALYSIS_SPECTRUMANALYZER_H

#include <string>

#include "DiaGraph.h"
#include "LevelScheme.h"

#include <TFile.h>
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
    TH2D gg;
    TGraph effgraph;
    //void SetupLevelSchemes();
    LevelScheme levelscheme;
    std::vector<std::pair<const Gamma*,const Gamma*>> gamma_gamma;
    const double fit_interval;
    const double proj_interval;
    unsigned fit_counter;
    double GetPeakIntegral(TH1D&, const double&);
};


#endif //EFFANALYSIS_SPECTRUMANALYZER_H
