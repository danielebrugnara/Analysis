#pragma once

#include <vector>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

class AgataProcessing
{
public:
    AgataProcessing();
    ~AgataProcessing();

    inline double CorrectDoppler(int, double);

private:
    unsigned long ref_ts;
    struct Gamma
    {
        const unsigned int multiplicity;
        bool in_coincidence;
        std::vector<double> E;
        std::vector<double> EDC;
        std::vector<TVector3> Pos;
        std::vector<TLorentzVector> Pgamma;
        Gamma(const unsigned int multiplicity)
            : multiplicity(multiplicity),
              in_coincidence(false)
        {

            E.resize(multiplicity);
            EDC.resize(multiplicity);
            Pos.resize(multiplicity);
            Pgamma.resize(multiplicity);
        };
    };

    const double z_shift;

public:
    inline void Process(){
        if (**data->AddTS - **data->LTS > ref_ts - 5 &&
            **data->AddTS - **data->LTS < ref_ts + 5)
            gammaray->in_coincidence = true;
        for (int ii = 0; ii < **(data->nbAdd); ++ii)
        {
            ComputeDoppler(ii);
        }
    };

    void ComputeDoppler(int);

    struct Data
    {
        TTreeReaderValue<int> *nbAdd;
        TTreeReaderValue<unsigned long long> *TSHit;
        TTreeReaderValue<unsigned long long> *AddTS;
        TTreeReaderValue<unsigned long long> *LTS;
        TTreeReaderArray<float> *AddE;
        TTreeReaderArray<float> *AddX;
        TTreeReaderArray<float> *AddY;
        TTreeReaderArray<float> *AddZ;
        TLorentzVector *p4;
        Data(TTreeReaderValue<int> *nbAdd,
             TTreeReaderValue<unsigned long long> *TSHit,
             TTreeReaderValue<unsigned long long> *AddTS,
             TTreeReaderValue<unsigned long long> *LTS,
             TTreeReaderArray<float> *AddE,
             TTreeReaderArray<float> *AddX,
             TTreeReaderArray<float> *AddY,
             TTreeReaderArray<float> *AddZ,
             TLorentzVector *p4) : nbAdd(nbAdd),
                                   TSHit(TSHit),
                                   AddTS(AddTS),
                                   LTS(LTS),
                                   AddE(AddE),
                                   AddX(AddX),
                                   AddY(AddY),
                                   AddZ(AddZ),
                                   p4(p4){};
    };

    inline void SetData(Data const *data)
    {
        if (this->data != nullptr)
            delete this->data;
        if (this->gammaray != nullptr)
            delete this->gammaray;
        this->data = data;
        gammaray = new Gamma(**(data->nbAdd));
    }

private:
    Gamma *gammaray;
    Data const *data;

public:
    inline bool In_Coincidence()
    {
        if (gammaray)
            return gammaray->in_coincidence;
        else
            return false;
    };
    inline unsigned int Get_Mult() { return gammaray->multiplicity; };
    inline double Get_E(const int &i) { return gammaray->E[i]; };
    inline double Get_EDC(const int &i) { return gammaray->EDC[i]; };
};