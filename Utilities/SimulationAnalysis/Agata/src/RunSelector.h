#pragma once

#include <iostream>
#include <string>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"
#include "TGraph2D.h"


#include "Fitter.h"
#include "Selector.h"

class RunSelector{
    public:
    RunSelector(TTree&, std::string );
    void Run(const std::vector<double>&, const std::string&);
    std::vector<Fitter> fits;
    int nevts;
    TTree& the_tree;
    std::string out_file_name;
private:
    void GetIntegral(TObject* obj, const std::vector<double>&);
    std::unordered_map<std::string, std::vector<double>> SearchCentroid(TObject* obj, const std::vector<double>&);
    int peak_search_cnt = 0;
};
