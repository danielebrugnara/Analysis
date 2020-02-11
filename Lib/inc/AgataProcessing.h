#ifndef __AGATA_PROCESSING_H__
#define __AGATA_PROCESSING_H__

#include <vector>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

class AgataProcessing{
    public:
        AgataProcessing(unsigned long, const double);
        ~AgataProcessing();

        inline double CorrectDoppler(int, double);
    private:
        unsigned long ref_ts;
        struct Gamma
        {
            const unsigned int multiplicity;
            std::vector<double> E;
            std::vector<double> EDC;
            std::vector<TVector3> Pos;
            std::vector<TLorentzVector> Pgamma;
            Gamma(const unsigned int multiplicity)
                    :   multiplicity(multiplicity){
                
                    E.resize(multiplicity);
                    EDC.resize(multiplicity);
                    Pos.resize(multiplicity);
                    Pgamma.resize(multiplicity);
                };
        };

        inline void Process();
        const double z_shift;
        inline void ComputeDoppler(int);
        
    public:
        struct Data
        {
            TTreeReaderValue<int> *            nbAdd;
            TTreeReaderValue<unsigned long> *  TSHit;
            TTreeReaderValue<unsigned long> *  AddTS;
            TTreeReaderArray<float> *          AddE;
            TTreeReaderArray<float> *          AddX;
            TTreeReaderArray<float> *          AddY;
            TTreeReaderArray<float> *          AddZ;
            TLorentzVector *                   p4;
        };
        inline void SetData(Data const*);
        
    private:
        Gamma*  gammaray;
        Data const*   data;
};

#endif