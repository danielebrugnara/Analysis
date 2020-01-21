
#ifndef __MUGASTIDENTIFICATION_H__
#define __MUGASTIDENTIFICATION_H__

#include "Identification.h"

class MugastIdentification: public Identification{
    public:
        MugastIdentification();
        ~MugastIdentification();

        bool Initialize();
        struct Data{
            //
            Data(){};
        };
        //cutsss: std::array<int, 3>
    private:
        std::unordered_map<int, std::unordered_map<int, double>> mass;
        const Double_t AMU_TO_MEV{931.4936148};
        Data const *data;

    public:
        
};

#endif