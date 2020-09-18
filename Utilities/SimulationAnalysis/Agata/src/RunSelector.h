#pragma once

#include <iostream>
#include <string>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"

class RunSelector{
    public:
    RunSelector(TTree&, const std::string&);
    int nevts;
};
