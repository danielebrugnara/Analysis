#pragma once

#include <vector>
#include <TVector3.h>
#include <TObject.h>

class MugastData: public TObject{
    public:
    const unsigned int multiplicity;            //Number of particles detected
    std::vector<TVector3> Pos;                  //3D positions
    std::vector<TVector3> EmissionDirection;    //Depends on the target position
    std::vector<TVector3> TelescopeNormal;      //
    std::vector<int> SI_X;                      //Intrinsic X
    std::vector<int> SI_Y;                      //Intrinsic Y
    std::vector<double> SI_E;                   //Energy deposition measured
    std::vector<double> SI_E2;                  //Second layer energy deposition
    std::vector<double> E;                      //After E-Loss corrections
    std::vector<double> E_CM;                   //After E-Loss corrections
    std::vector<double> Ex;                     //Ex computed with reaction
    std::vector<double> E2;                     //Second layer
    std::vector<double> SI_T;                   //Un calibrated time
    std::vector<double> T;                      //Calibrated time
    std::vector<double> T2;                     //
    std::vector<double> MG;                     //
    std::vector<int> M;                         //Mass number
    std::vector<int> Z;                         //Z number
    std::vector<bool> Indentified;              //
    std::vector<std::string> Particle;          //

    MugastData(): multiplicity(0){};
    explicit MugastData(const unsigned int multiplicity)
            : multiplicity(multiplicity)
    {
        Pos.resize(multiplicity);
        EmissionDirection.resize(multiplicity);
        TelescopeNormal.resize(multiplicity);
        SI_X.resize(multiplicity);
        SI_Y.resize(multiplicity);
        SI_E.resize(multiplicity);
        SI_E2.resize(multiplicity);
        E.resize(multiplicity);
        E_CM.resize(multiplicity);
        Ex.resize(multiplicity);
        E2.resize(multiplicity);
        SI_T.resize(multiplicity);
        T.resize(multiplicity);
        T2.resize(multiplicity);
        MG.resize(multiplicity);
        M.resize(multiplicity);
        Z.resize(multiplicity);
        Indentified.resize(multiplicity);
        Particle.resize(multiplicity);
    };
    ClassDef(MugastData, 1);
}; 
