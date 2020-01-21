#ifndef __IDENTIFICATION_H__
#define __IDENTIFICATION_H__

#include <unordered_map>
#include <fstream>
#include <iostream>

#include <TCutG.h>
#include <TFile.h>
#include <TKey.h>

class Identification {
   public:
    Identification();
    ~Identification();

    void LoadCuts(std::string);
   private:

   protected:
    //TODO: use templates to include SetData here
    std::unordered_map<std::string, TCutG*> cuts;
    std::unordered_map<std::string, std::unordered_map<std::string, TCutG*>> cut_type;
};

#endif