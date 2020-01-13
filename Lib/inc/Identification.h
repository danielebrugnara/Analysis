#ifndef __IDENTIFICATION_H__
#define __IDENTIFICATION_H__

#include <unordered_map>
#include <fstream>

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
    std::unordered_map<std::string, TCutG*> cuts;
    std::unordered_map<std::string, std::unordered_map<std::string, TCutG*>> cut_type;
};

#endif