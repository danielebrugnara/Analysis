#pragma once

#include <TVector3.h>
#include <TLorentzVector.h>

class AgataData: public TObject{
    public:
    const unsigned int multiplicity;
    unsigned long TimeDiff;
    bool in_coincidence;
    std::vector<double> E;
    std::vector<double> EDC;
    std::vector<double> EDCMidTarget;
    std::vector<TVector3> Pos;
    std::vector<TLorentzVector> Pgamma;
    std::vector<TLorentzVector> PgammaMidTarget;
    AgataData(): multiplicity(0), in_coincidence(false){};
    explicit AgataData(const unsigned int multiplicity)
            : multiplicity(multiplicity),
              TimeDiff(0),
              in_coincidence(false)
    {

        E.resize(multiplicity);
        EDC.resize(multiplicity);
        EDCMidTarget.resize(multiplicity);
        Pos.resize(multiplicity);
        Pgamma.resize(multiplicity);
        PgammaMidTarget.resize(multiplicity);
    };
    ClassDef(AgataData, 3);
};