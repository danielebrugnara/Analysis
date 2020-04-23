
#ifndef __VAMOSIDENTIFICATION_H__
#define __VAMOSIDENTIFICATION_H__

#include <TLorentzVector.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "Identification.h"
#include "Interpolation.h"

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
    struct Fragment
    {
        double En;         //Energy IC
        double D_En;       //DE 1+2
        double D_En2;      //DE 1
        double Path;       //Reconstructed path
        double T;          //Time
        double V;          //Velocity
        double Beta;       //Beta
        double Gamma;      //Gamma
        double M_Q;        //Mass over charge
        double M;          //Mass
        double Charge;     //Charge state
        bool Identified;   //Positive identification
        TLorentzVector p4; //4 momentum of recoil
        int id_M;
        int id_Z;
        int id_Q;
        Fragment() : En(0),
                     D_En(0),
                     D_En2(0),
                     Path(0),
                     T(0),
                     V(0),
                     Beta(0),
                     Gamma(0),
                     M_Q(0),
                     M(0),
                     Charge(0),
                     Identified(false),
                     p4(),
                     id_M(0),
                     id_Z(0),
                     id_Q(0){};
    };
    Data const *data;
    Fragment *fragment;

    Interpolation *FP_time_interpolation;

    const double IC_threashold{0.1};

public:
    //Various computed fragment properties
    inline double Get_En() { return fragment->En; };
    inline double Get_D_En() { return fragment->D_En; };
    inline double Get_D_En2() { return fragment->D_En2; };
    inline double Get_Path() { return fragment->Path; };
    inline double Get_T() { return fragment->T; };
    inline double Get_V() { return fragment->V; };
    inline double Get_Beta() { return fragment->Beta; };
    inline double Get_Gamma() { return fragment->Gamma; };
    inline double Get_M_Q() { return fragment->M_Q; };
    inline double Get_M() { return fragment->M; };
    inline double Get_Charge() { return fragment->Charge; };
    inline double Identified() { return fragment->Identified; };
    inline TLorentzVector *Get_p4() { return &(fragment->p4); };
    inline unsigned int Get_id_Z() { return fragment->id_Z; };
    inline unsigned int Get_id_M() { return fragment->id_M; };
    inline unsigned int Get_id_Q() { return fragment->id_Q; };

    //////////////////////////////////////////////////////////////////
    //Inline Functions implementation (required to be in header file)
public:
    inline void SetData(Data const *data)
    {
#ifdef VERBOSE_DEBUG
        std::cout << "------------>VamosIdentification::SetData()\n";
#endif
        //TODO: will this cause memory leaks?
        if (this->data != nullptr)
            delete this->data;
        if (this->fragment != nullptr)
            delete this->fragment;
#ifdef VERBOSE_DEBUG
        std::cout << "------------>finished : pointers deleted\n";
#endif
        fragment = new Fragment();
        this->data = data;
    }

    bool Identify();

private:
    double GetShift();

    inline double GetFPTime()
    {
        return 540.5 * (**data->AGAVA_VAMOSTS < 104753375647998) + 537.9 * (**data->AGAVA_VAMOSTS >= 104753375647998) - 2. * **data->T_FPMW_CATS2_C + GetShift();
    }
};

#endif
