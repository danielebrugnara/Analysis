#ifndef RUN_SELECTOR_H
#define RUN_SELECTOR_H

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

#include "Selector.h"

class RunSelector{
    public:
    RunSelector(std::string);
    int nevts;
};
#endif
