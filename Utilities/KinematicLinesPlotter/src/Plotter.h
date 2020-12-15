#pragma once

#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

#include <TStyle.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TCanvas.h>

class Plotter{
public:
    Plotter()=default;
    void plotLines();
};