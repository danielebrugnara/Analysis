#pragma once

#include <iostream>
#include <string>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"


#include "Fitter.h"
#include "Selector.h"

class RunSelector{
    public:
    RunSelector(TTree&, std::string );
    void Run(const std::vector<double>&);
    std::vector<Fitter> fits;
    int nevts;
    TTree& the_tree;
    std::string out_file_name;
};
