#pragma once

#include <vector>
#include <TVector3.h>
#include <TObject.h>
#include <MugastData.h>

class Must2Data: public MugastData{
    public:
    std::vector<double> CsI_E;                   //Energy deposition in scintillator measured
    std::vector<double> CsI_T;                   //Time in scintillator measured

    Must2Data() = default;

    explicit Must2Data(const unsigned int multiplicity)
            :   MugastData(multiplicity)
    {
        CsI_E.resize(multiplicity);
        CsI_T.resize(multiplicity);
    };
    ClassDef(Must2Data, 2);
}; 
