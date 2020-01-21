
#ifndef __MUGASTIDENTIFICATION_H__
#define __MUGASTIDENTIFICATION_H__

#include "Identification.h"
#include "Interpolation.h"

#include <TTreeReaderValue.h>

#include "TMugastPhysics.h"

#include <array>

class MugastIdentification: public Identification{
    public:
        MugastIdentification();
        ~MugastIdentification();
        bool Initialize();

        struct Data{
            TTreeReaderValue<TMugastPhysics> *Mugast;
            Data(TTreeReaderValue<TMugastPhysics> *Mugast):
                    Mugast(Mugast){};
        };
        std::array<int, 6> cuts_MG;
        std::array<int, 3> cuts_M;
        std::array<int, 2> cuts_Z;
    private:
        struct Fragment{
            const unsigned int multiplicity;
            std::vector<double> X;
            std::vector<double> Y;
            std::vector<double> Z;
            std::vector<double> E;
            std::vector<double> T;
            std::vector<double> MG;
            Fragment(const unsigned int multiplicity) : multiplicity(multiplicity){
                X.resize(multiplicity);
                Y.resize(multiplicity);
                Z.resize(multiplicity);
                E.resize(multiplicity);
                T.resize(multiplicity);
                MG.resize(multiplicity);
            };
        };
        Data const *data;
        Fragment *fragment;

        Interpolation *gas_thickness;
        Interpolation *havar_angle;
        std::unordered_map<int, std::unordered_map<int, double>> mass;
        const Double_t AMU_TO_MEV{931.4936148};

    public:
        inline void SetData(Data const *data){
            if (this->data != nullptr) delete this->data;
            if (this->fragment != nullptr) delete this->fragment;
            this->data = data;
            fragment = new Fragment((**(data->Mugast)).DSSD_E.size());
        }
        inline bool Identify(){
            //for (const auto &particle_search : )
        }
};

#endif