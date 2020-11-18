#pragma once

#include <TVector3.h>
#include <TLorentzVector.h>

class AgataData: public TObject{
    public:
    const unsigned int multiplicity;
    bool in_coincidence;
    std::vector<double> E;
    std::vector<double> EDC;
    std::vector<TVector3> Pos;
    std::vector<TLorentzVector> Pgamma;
    AgataData(): multiplicity(0), in_coincidence(false){};
    explicit AgataData(const unsigned int multiplicity)
            : multiplicity(multiplicity),
              in_coincidence(false)
    {

        E.resize(multiplicity);
        EDC.resize(multiplicity);
        Pos.resize(multiplicity);
        Pgamma.resize(multiplicity);
    };
    ClassDef(AgataData, 1);
};