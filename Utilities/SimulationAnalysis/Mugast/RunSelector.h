#ifndef RUN_SELECTOR_H
#define RUN_SELECTOR_H

#include <iostream>
#include <string>
#include <map>

#include "TFile.h"
#include "TTree.h"

#include "Selector.h"
#include "SelectorData.h"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>



#include <fstream>

class RunSelector{
    public:
    RunSelector(std::string, std::string);
};
#endif
