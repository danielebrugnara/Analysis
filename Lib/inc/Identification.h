#ifndef __IDENTIFICATION_H__
#define __IDENTIFICATION_H__

#include <map>

#include <TCutG.h>

class Identification {
   public:
    Identification();
    ~Identification();
   private:
    std::map<std::string, TCutG*> cuts;
};

#endif