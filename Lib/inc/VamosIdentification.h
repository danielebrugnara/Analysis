#pragma once

#include <TLorentzVector.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "Identification.h"
#include "Interpolation.h"
#include "VamosData.h"
#include "ReactionFragment.h"

#include <array>
#include <sstream>

class VamosIdentification : public Identification
{
public:
    VamosIdentification();
    ~VamosIdentification();
    bool initialize();

    struct Data
    {
        TTreeReaderArray<float> *IC;
        TTreeReaderValue<float> *Path;
        TTreeReaderValue<float> *Brho;
        TTreeReaderValue<float> *Xf;
        TTreeReaderValue<float> *Yf;
        TTreeReaderValue<float> *Pf;
        TTreeReaderValue<float> *Tf;
        TTreeReaderValue<float> *ThetaL;
        TTreeReaderValue<float> *PhiL;
        TTreeReaderValue<long long unsigned int> *AGAVA_VAMOSTS;
        TTreeReaderValue<float> *T_FPMW_CATS2_C;
        Data(TTreeReaderArray<float> *IC,
             TTreeReaderValue<float> *Path,
             TTreeReaderValue<float> *Brho,
             TTreeReaderValue<float> *Xf,
             TTreeReaderValue<float> *Yf,
             TTreeReaderValue<float> *Pf,
             TTreeReaderValue<float> *Tf,
             TTreeReaderValue<float> *ThetaL,
             TTreeReaderValue<float> *PhiL,
             TTreeReaderValue<long long unsigned int> *AGAVA_VAMOSTS,
             TTreeReaderValue<float> *T_FPMW_CATS2_C) : IC(IC),
                                                        Path(Path),
                                                        Brho(Brho),
                                                        Xf(Xf),
                                                        Yf(Yf),
                                                        Pf(Pf),
                                                        Tf(Tf),
                                                        ThetaL(ThetaL),
                                                        PhiL(PhiL),
                                                        AGAVA_VAMOSTS(AGAVA_VAMOSTS),
                                                        T_FPMW_CATS2_C(T_FPMW_CATS2_C){};
    };
    std::array<int, 3> cutsZ;
    //std::array<int, 5> cutsM;
    //std::array<int, 9> cuts_Q;
    std::array<int, 4> cutsM;
    //std::array<int, 8> cuts_Q;

private:
    std::unordered_map<int, std::unordered_map<int, double>> mass;
    //const Double_t AMU_TO_MEV{931.4936148};
    void readFpTimeShifts();
    std::vector<std::pair<double, double>> timeShifts; //Xf-max, shift
    Data const *data;
    VamosData fragment;

    Interpolation *fpTimeInterpolation;

    const double icThreashold{0.1};

public:
    //Various computed fragment properties

    inline double getEn()            const { return fragment.En; };
    inline double getDEn()          const { return fragment.D_En; };
    inline double getDEn2()         const { return fragment.D_En2; };
    inline double getPath()          const { return fragment.Path; };
    inline double getT()             const { return fragment.T; };
    inline double getV()             const { return fragment.V; };
    inline double getBeta()          const { return fragment.Beta; };
    inline double getGamma()         const { return fragment.Gamma; };
    inline double getMQ()           const { return fragment.M_Q; };
    inline double getM()             const { return fragment.M; };
    inline double getCharge()        const { return fragment.Charge; };
    inline double identified()        const { return fragment.Identified; };
    inline TLorentzVector *getP4()         { return &(fragment.p4); };
    inline TLorentzVector *getP4MidTarget()         { return &(fragment.p4MidTarget); };
    inline unsigned int getIdZ()    const { return fragment.id2_Z; };
    inline unsigned int getIdM()    const { return fragment.id_M; };
    inline unsigned int getIdQ()    const { return fragment.id_Q; };

    inline void setP4MidTarget(const double& beamEnergy) {
        fragment.p4MidTarget = *getP4();
        if (fragment.p4MidTarget.M() > 0) {
            TVector3 p = fragment.p4MidTarget.Vect();
            double Etot = beamEnergy + fragment.p4MidTarget.M();
            p.SetMag(sqrt(pow(Etot, 2) - fragment.p4MidTarget.M2()));
            fragment.p4MidTarget.SetVectM(p, fragment.p4MidTarget.M());
        }
    };

    inline VamosData& getData()        {return fragment;};

    double getEnFromBrho();
    //////////////////////////////////////////////////////////////////
    //Inline Functions implementation (required to be in header file)
public:
    inline void SetData(Data const *data)
    {
        delete this->data;
        this->data = data;

        //The following constructs at the same address
        fragment.~VamosData();
        new(&fragment) VamosData();
    }

    bool identify();

private:
    double getShift();

    inline double getFpTime()
    {
        return 540.5 * (**data->AGAVA_VAMOSTS < 104753375647998) + 537.9 * (**data->AGAVA_VAMOSTS >= 104753375647998) - 2. * **data->T_FPMW_CATS2_C +
                getShift();
    }
};