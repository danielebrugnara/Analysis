#pragma once

#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

#include <TFile.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TCanvas.h>

class Plotter{
public:
    Plotter();
    void plotVamosAcceptance(const std::string&);
    void plotMugastAcceptance(const std::string&);
private:
    TStyle style;
    std::map<std::string, ReactionReconstruction2body<long double>*> reactions;
    double beam_energy;
    //std::vector<int> colors{46, 36, 32, 41, 226, 221};
    std::vector<int> colors{1, 2, 3, 4, 6, 7, 9};
};