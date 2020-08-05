#ifndef RUN_SELECTOR_H
#define RUN_SELECTOR_H

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TVector2.h"

#include "Selector.h"

class RunSelector{
public:
    RunSelector(std::string);
    std::string GetFileName() const;
private:
    std::string file_name;
};
#endif
