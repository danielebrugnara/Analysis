#ifndef __IDENTIFICATION_H__
#define __IDENTIFICATION_H__

#include <map>

#include <TCutG.h>
#include <TFile.h>
#include <TKey.h>

class Identification {
   public:
    Identification();
    ~Identification();

    void LoadCuts(std::string);
   private:
    std::map<std::string, TCutG*> cuts;
};

#endif