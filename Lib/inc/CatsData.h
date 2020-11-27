#pragma once

#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

class CatsData: public TObject{
public:
    const unsigned int multiplicity;            //Number of particles detected
    std::vector<TVector3> Pos;
    //std::vector<TLorentzVector> Beam;
    std::vector<TVector2> Charge;
    CatsData(): multiplicity(0){};

    explicit CatsData(const unsigned int multiplicity)
                    :multiplicity(multiplicity){
        Pos.resize(multiplicity);
        //Beam.resize(multiplicity);
        Charge.resize(multiplicity);
    }
    ClassDef(CatsData, 1);
}; 
