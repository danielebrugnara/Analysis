#pragma once

#include <vector>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <AgataData.h>

class AgataProcessing
{
public:
    AgataProcessing();
    ~AgataProcessing();

    inline double CorrectDoppler(int, double);

private:
    unsigned long   ref_ts;
    const double    z_shift;

public:
    inline void Process(){
        if (**data->AddTS - **data->LTS > ref_ts - 5 &&
            **data->AddTS - **data->LTS < ref_ts + 5)
            gammaray.in_coincidence = true;
        for (int ii = 0; ii < **(data->nbAdd); ++ii){
            ComputeDoppler(ii);
        }
    };

    void ComputeDoppler(int);

    struct Data{
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

    inline void SetData(Data const *data){
        delete this->data;
        this->data = data;

        //The following constructs at the same address
        gammaray.~AgataData();
        new (&gammaray) AgataData(**(data->nbAdd));
    }

private:
    AgataData gammaray;
    Data const *data;

public:
    inline bool In_Coincidence()            const { return gammaray.in_coincidence; };
    inline unsigned int Get_Mult()          const { return gammaray.multiplicity; };
    inline double Get_E(const int &i)       const { return gammaray.E[i]; };
    inline double Get_EDC(const int &i)     const { return gammaray.EDC[i]; };
    inline AgataData& Get_Data()                  { return gammaray; };
};