#pragma once

#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

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
    void plotVamosAcceptance();
    void plotMugastAcceptance();
private:
    TStyle style;
    std::map<std::string, ReactionReconstruction2body<long double>*> reactions;
    double beam_energy;
};