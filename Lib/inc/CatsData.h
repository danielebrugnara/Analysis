#pragma once

#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

class CatsData: public TObject{
public:
    const unsigned int multiplicity;            //Number of particles detected
    std::vector<TVector3> Pos;
    std::vector<TVector3> PosOnTarget;
    std::vector<TVector3> Dir;
    //std::vector<TLorentzVector> Beam;
    std::vector<TVector3> Charge;
    CatsData(): multiplicity(0){};

    explicit CatsData(const unsigned int multiplicity, const unsigned int options)
                    :multiplicity(multiplicity==0 ? 1 : multiplicity){
        Pos.resize(this->multiplicity);
        PosOnTarget.resize(options);
        Dir.resize(options);
        //Beam.resize(multiplicity);
        Charge.resize(this->multiplicity);
    }
    ClassDef(CatsData, 4);
}; 
