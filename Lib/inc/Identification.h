#pragma once

#include <unordered_map>
#include <fstream>
#include <iostream>

#include <TCutG.h>
#include <TFile.h>
#include <TKey.h>

class Identification
{
public:
    Identification();
    ~Identification();

    void LoadCuts(const std::string&);

private:
protected:
    //TODO: use templates to include setData here
    std::unordered_map<std::string, TCutG *> cuts;
    std::unordered_map<std::string, std::unordered_map<std::string, TCutG *>> cut_type;
};