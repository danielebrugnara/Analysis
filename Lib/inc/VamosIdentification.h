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
    bool Initialize();

    struct Data
    {
        TTreeReaderArray<float> *IC;
        TTreeReaderValue<float> *Path;
        TTreeReaderValue<float> *Brho;
        TTreeReaderValue<float> *Xf;
        TTreeReaderValue<float> *ThetaL;
        TTreeReaderValue<float> *PhiL;
        TTreeReaderValue<long long unsigned int> *AGAVA_VAMOSTS;
        TTreeReaderValue<float> *T_FPMW_CATS2_C;
        Data(TTreeReaderArray<float> *IC,
             TTreeReaderValue<float> *Path,
             TTreeReaderValue<float> *Brho,
             TTreeReaderValue<float> *Xf,
             TTreeReaderValue<float> *ThetaL,
             TTreeReaderValue<float> *PhiL,
             TTreeReaderValue<long long unsigned int> *AGAVA_VAMOSTS,
             TTreeReaderValue<float> *T_FPMW_CATS2_C) : IC(IC),
                                                        Path(Path),
                                                        Brho(Brho),
                                                        Xf(Xf),
                                                        ThetaL(ThetaL),
                                                        PhiL(PhiL),
                                                        AGAVA_VAMOSTS(AGAVA_VAMOSTS),
                                                        T_FPMW_CATS2_C(T_FPMW_CATS2_C){};
    };
    std::array<int, 3> cuts_Z;
    //std::array<int, 5> cuts_M;
    //std::array<int, 9> cuts_Q;
    std::array<int, 4> cuts_M;
    std::array<int, 8> cuts_Q;

private:
    std::unordered_map<int, std::unordered_map<int, double>> mass;
    const Double_t AMU_TO_MEV{931.4936148};
    void ReadFPTimeShifts();
    std::vector<std::pair<double, double>> TimeShifts; //Xf-max, shift
    Data const *data;
    VamosData fragment;

    Interpolation *FP_time_interpolation;

    const double IC_threashold{0.1};

public:
    //Various computed fragment properties

    inline double Get_En()            const { return fragment.En; };
    inline double Get_D_En()          const { return fragment.D_En; };
    inline double Get_D_En2()         const { return fragment.D_En2; };
    inline double Get_Path()          const { return fragment.Path; };
    inline double Get_T()             const { return fragment.T; };
    inline double Get_V()             const { return fragment.V; };
    inline double Get_Beta()          const { return fragment.Beta; };
    inline double Get_Gamma()         const { return fragment.Gamma; };
    inline double Get_M_Q()           const { return fragment.M_Q; };
    inline double Get_M()             const { return fragment.M; };
    inline double Get_Charge()        const { return fragment.Charge; };
    inline double Identified()        const { return fragment.Identified; };
    inline TLorentzVector *Get_p4()         { return &(fragment.p4); };
    inline unsigned int Get_id_Z()    const { return fragment.id_Z; };
    inline unsigned int Get_id_M()    const { return fragment.id_M; };
    inline unsigned int Get_id_Q()    const { return fragment.id_Q; };

    inline VamosData& Get_Data()        {return fragment;};

    double Get_EnFromBrho();
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

    bool Identify();

private:
    double GetShift();

    inline double GetFPTime()
    {
        return 540.5 * (**data->AGAVA_VAMOSTS < 104753375647998) + 537.9 * (**data->AGAVA_VAMOSTS >= 104753375647998) - 2. * **data->T_FPMW_CATS2_C + GetShift();
    }
};