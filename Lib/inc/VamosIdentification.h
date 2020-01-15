
#ifndef __VAMOSIDENTIFICATION_H__
#define __VAMOSIDENTIFICATION_H__

#include <TLorentzVector.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include "Identification.h"
#include "Interpolation.h"

#include <array>
#include <sstream>

class VamosIdentification : public Identification {
   public:
    VamosIdentification();
    ~VamosIdentification();
    bool Initialize();

    struct Data {
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
    std::array<int, 5> cuts_M;
    std::array<int, 9> cuts_Q;

    std::unordered_map<int, std::unordered_map<int, double>> mass;

   private:
    const Double_t AMU_TO_MEV{931.4936148};
    void ReadFPTimeShifts();
    std::vector<std::pair<double, double>> TimeShifts;  //Xf-max, shift
    struct Fragment {
        double En;          //Energy IC
        double D_En;        //DE 1+2
        double D_En2;       //DE 1
        double Path;        //Reconstructed path
        double T;           //Time
        double V;           //Velocity
        double Beta;        //Beta
        double Gamma;       //Gamma
        double M_Q;         //Mass over charge
        double M;           //Mass
        double Charge;      //Charge state
        bool Identified;    //Positive identification
        TLorentzVector p4;  //4 momentum of recoil
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
    inline double Get_id_Z() { return fragment->id_Z; };
    inline double Get_id_M() { return fragment->id_M; };
    inline double Get_id_Q() { return fragment->id_Q; };

    //////////////////////////////////////////////////////////////////
    //Inline Functions implementation (required to be in header file)
   public:
    inline void SetData(Data const *data) {
#ifdef VERBOSE_DEBUG
        std::cout << "------------>VamosIdentification::SetData()\n";
#endif
        //TODO: will this cause memory leaks?
        //if (this->data)     delete this->data;
        //if (this->fragment!=nullptr) delete (this->fragment);
#ifdef VERBOSE_DEBUG
        std::cout << "------------>finished : pointers deleted\n";
#endif
        fragment = new Fragment();
        this->data = data;
    }
    inline bool Identify() {
        fragment->En        = ((*data->IC)[0] > IC_threashold) * ((*data->IC)[0] + ((*data->IC)[1] > IC_threashold) * ((*data->IC)[1] + ((*data->IC)[2] > IC_threashold) * ((*data->IC)[2] + ((*data->IC)[3] > IC_threashold) * ((*data->IC)[3] + ((*data->IC)[4] > IC_threashold) * ((*data->IC)[4] + ((*data->IC)[5] > IC_threashold) * ((*data->IC)[5]))))));
        fragment->D_En      = ((*data->IC)[0] > IC_threashold) * ((*data->IC)[0] + ((*data->IC)[1] > IC_threashold) * ((*data->IC)[1]));
        fragment->D_En2     = (*data->IC)[0] * ((*data->IC)[1] > IC_threashold);

        //fragment->T = 540.5 * (data->AGAVA_VAMOSTS < 104753375647998) + 537.9 * (data->AGAVA_VAMOSTS >= 104753375647998) - 2. * (data->T_FPMW_CATS2_C) + 2.7 * (MW_N[0] == 16) + 2.7 * (MW_N[0] == 15) + 2.9 * (MW_N[0] == 14) + 2.9 * (MW_N[0] == 13) + 2.4 * (MW_N[0] == 12) + 1.3 * (MW_N[0] == 11) + 1.5 * (MW_N[0] == 10) + 1.6 * (MW_N[0] == 9) - 0.6 * (MW_N[0] == 8) + 2.5 * (MW_N[0] == 7) + 2. * (MW_N[0] == 6) + 1.6 * (MW_N[0] == 5) + 1.1 * (MW_N[0] == 4) - 0.6 * (MW_N[0] == 3) - 1.2 * (MW_N[0] == 2) - 4.0 * (MW_N[0] == 1);
        fragment->T         = GetFPTime();
        fragment->Path      = **data->Path + 5;

        //Computing the basic identifiaction
        fragment->V         = fragment->Path / fragment->T;
        fragment->Beta      = fragment->V / 29.9792;
        fragment->Gamma     = 1. / sqrt(1.0 - fragment->Beta * fragment->Beta);
        fragment->M         = (fragment->En) / 931.5016 / (fragment->Gamma - 1.);
        //mM2 = 18./20.8*(mE2)/931.5016/(mGamma2-1.);
        fragment->M_Q       = **data->Brho / 3.105 / fragment->Beta / fragment->Gamma;
        fragment->Charge    = fragment->M / fragment->M_Q;

        //dE - E identification
        for (const auto &z_search : cut_type.at("dE2_E")) {
            if (z_search.second->IsInside(fragment->En, fragment->D_En2)) {
                //Z format dE2_E_Z18
                if (fragment->id_Z == 0)
                    fragment->id_Z = stoi(
                        z_search.first.substr(z_search.first.find_last_of("_Z") + 1));
                else
                    throw std::runtime_error("Overlapping Z gates\n");
            }
        }
        if (fragment->id_Z == 0) return false;

        //MQ - Q identification
        for (const auto &mq_search : cut_type.at("MQ_Q")) {
            if (mq_search.second->IsInside(fragment->M_Q, fragment->Charge)) {
                if (fragment->id_M == 0 && fragment->id_Q == 0) {
                    fragment->id_M = stoi(mq_search.first.substr(mq_search.first.find_last_of("M") + 1, 2));
                    fragment->id_Q = stoi(mq_search.first.substr(mq_search.first.find_last_of("_") + 2));
                } else
                    throw std::runtime_error("Overlapping M_Q gates");
            }
        }
        if (fragment->id_M == 0 || fragment->id_Q == 0) return false;

        //Lorentzvector computation
        fragment->p4.SetT(mass[fragment->id_M][fragment->id_Z]);
        TVector3 v4(0, 0, fragment->Beta);
        v4.SetMagThetaPhi(fragment->Beta, **data->ThetaL, **data->PhiL);
        fragment->p4.Boost(v4);

        return fragment->Identified = true;
    }

   private:
    inline double GetShift() {
        double shift{0};
        for (const auto &it : TimeShifts) {
            if (**data->Xf < it.first)
                shift = it.first;
            else
                break;
        }
        return shift;
    }

    inline double GetFPTime() {
        //    return 540.5*(data->AGAVA_VAMOSTS<104753375647998)+537.9*(data->AGAVA_VAMOSTS>=104753375647998) -2.*data->T_FPMW_CATS2_C + GetShift();
        return 0;
    }
};

#endif
