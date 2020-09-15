#pragma once

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"

class RunSelector{
    public:
    RunSelector(std::string);
    int nevts;
    std::string file_out_name;
};
