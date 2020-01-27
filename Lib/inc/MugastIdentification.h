
#ifndef __MUGASTIDENTIFICATION_H__
#define __MUGASTIDENTIFICATION_H__

#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <array>

#include "Identification.h"
#include "Interpolation.h"
#include "Calibration.h"
#include "TMugastPhysics.h"

//TODO: generate internal library for the following headers
#include "NPEnergyLoss.h"

class MugastIdentification : public Identification {
   public:
    MugastIdentification();
    ~MugastIdentification();
    bool Initialize();

    static constexpr int n_detectors = 6;
    static constexpr int n_strips = 128;

    struct Data {
        TTreeReaderValue<TMugastPhysics> *Mugast;
        Data(TTreeReaderValue<TMugastPhysics> *Mugast) : Mugast(Mugast){};
    };
    std::array<int, n_detectors> cuts_MG;
    std::array<int, 3> cuts_M;
    std::array<int, 2> cuts_Z;
    std::array<std::string, 3> particles;
    std::array<std::string, 2> strips;

   private:
    struct Fragment {
        const unsigned int multiplicity;
        std::vector<TVector3> Pos;
        std::vector<int> SI_X;
        std::vector<int> SI_Y;
        std::vector<double> SI_E;
        std::vector<double> SI_E2;
        std::vector<double> E;
        std::vector<double> E2;
        std::vector<double> T;
        std::vector<double> T2;
        std::vector<double> MG;
        std::vector<int> M;
        std::vector<int> Z;
        std::vector<double> Ex;
        Fragment(const unsigned int multiplicity) : multiplicity(multiplicity) {
            Pos.resize(multiplicity);
            SI_X.resize(multiplicity);
            SI_Y.resize(multiplicity);
            SI_E.resize(multiplicity);
            SI_E2.resize(multiplicity);
            E.resize(multiplicity);
            E2.resize(multiplicity);
            T.resize(multiplicity);
            T2.resize(multiplicity);
            MG.resize(multiplicity);
            M.resize(multiplicity);
            Z.resize(multiplicity);
        };
    };
    Data const *data;
    Fragment *fragment;

    //Interpolations
    Interpolation *gas_thickness;
    Interpolation *havar_angle;

    //Calibrations
    std::unordered_map<int, Calibration *> calibrations_TY;

    //Physics
    const Double_t AMU_TO_MEV{931.4936148};
    std::unordered_map<int, std::unordered_map<int, double>> mass;

    //Energy Loss
    std::unordered_map<std::string,std::map<std::string, NPL::EnergyLoss *>> energy_loss;

    bool with_cuts;

    bool InitializeCuts();
    bool InitializeCalibration();
    bool InitializeELoss();

   public:
    inline void SetData(Data const *data) {
        if (this->data != nullptr)
            delete this->data;
        if (this->fragment != nullptr)
            delete this->fragment;
        this->data = data;
        fragment = new Fragment((**(data->Mugast)).DSSD_E.size());
    };

    inline bool Identify() {
#ifdef VERBOSE_DEBUG
        std::cout << "------------>MugastIdentification::Identify()\n";
#endif
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            fragment->Pos[ii] = TVector3((**(data->Mugast)).PosX[ii],
                                         (**(data->Mugast)).PosY[ii],
                                         (**(data->Mugast)).PosZ[ii]);
            fragment->SI_E[ii] = (**(data->Mugast)).DSSD_E[ii];
            fragment->SI_E2[ii] = (**(data->Mugast)).SecondLayer_E[ii];
            fragment->SI_X[ii] = (**(data->Mugast)).DSSD_X[ii];
            fragment->SI_Y[ii] = (**(data->Mugast)).DSSD_Y[ii];
            fragment->T[ii] = (**(data->Mugast)).DSSD_T[ii];
            fragment->T2[ii] = (**(data->Mugast)).SecondLayer_T[ii];
            fragment->MG[ii] = (**(data->Mugast)).TelescopeNumber[ii];
        }
        try {
            for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
                bool already_id = false;
                for (const auto &cut_it : particles) {
                    if (with_cuts && cut_type["E_TOF"].at("E_TOF_" + cut_it + "_MG" + std::to_string(fragment->MG[ii]))->IsInside(fragment->SI_E[ii], fragment->T[ii])) {
                        if (already_id) throw std::runtime_error("Overlapping MUGAST E TOF gates");

                        fragment->M[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("m") + 1,
                                                                  cut_it.find_first_of("_") - cut_it.find_first_of("m") - 1));

                        fragment->Z[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("z") + 1,
                                                                  cut_it.find_first_of("_") - cut_it.find_first_of("z") - 1));

                        already_id = true;
                    }
                }
                if (!already_id) {
                    fragment->M[ii] = 0;
                    fragment->Z[ii] = 0;
                }
            };
        } catch (const std::out_of_range &err) {
            std::cerr << "Mugast cuts not found\n";
            with_cuts = false;
        }
        return true;
        //TODO: reconstruct E here

        //TODO: compute Ex here
    }

    inline int          Get_Mult()              { return fragment->multiplicity; };
    inline TVector3 *   Get_Pos(const int &i)   { return &(fragment->Pos[i]); };
    inline int          Get_SI_X(const int &i)  { return fragment->SI_X[i]; };
    inline int          Get_SI_Y(const int &i)  { return fragment->SI_Y[i]; };
    inline double       Get_E(const int &i)     { return fragment->E[i]; };
    inline double       Get_SI_E(const int &i)  { return fragment->SI_E[i]; };
    inline double       Get_T2(const int &i)    { return fragment->T2[i]; };
    inline double       Get_MG(const int &i)    { return fragment->MG[i]; };
    inline double       Get_M(const int &i)     { return fragment->M[i]; };
    inline double       Get_Z(const int &i)     { return fragment->Z[i]; };

    inline double       Get_T(const int &i)     { if (calibrations_TY[fragment->MG[i]]==nullptr) 
                                                    return fragment->T[i]; 
                                                  else
                                                    return calibrations_TY[fragment->MG[i]]
                                                            ->Evaluate(fragment->T[i], fragment->SI_Y[i]);};
};

#endif