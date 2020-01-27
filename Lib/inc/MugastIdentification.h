
#ifndef __MUGASTIDENTIFICATION_H__
#define __MUGASTIDENTIFICATION_H__

#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <array>

#include "Calibration.h"
#include "Identification.h"
#include "Interpolation.h"
#include "TMugastPhysics.h"

//TODO: generate internal library for the following headers
#include "NPEnergyLoss.h"

class MugastIdentification : public Identification {
   public:
    MugastIdentification();
    ~MugastIdentification();
    bool Initialize(const double &, const TVector3 &);  //Beam energy in MeV

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
    double beam_energy;
    TVector3 target_pos;
    struct Fragment {
        const unsigned int multiplicity;
        std::vector<TVector3> Pos;
        std::vector<int> SI_X;
        std::vector<int> SI_Y;
        std::vector<double> SI_E;
        std::vector<double> SI_E2;
        std::vector<double> E;
        std::vector<double> E2;
        std::vector<double> SI_T;
        std::vector<double> T;
        std::vector<double> T2;
        std::vector<double> MG;
        std::vector<int> M;
        std::vector<int> Z;
        std::vector<double> Ex;
        std::vector<bool> Indentified;
        Fragment(const unsigned int multiplicity) : multiplicity(multiplicity) {
            Pos.resize(multiplicity);
            SI_X.resize(multiplicity);
            SI_Y.resize(multiplicity);
            SI_E.resize(multiplicity);
            SI_E2.resize(multiplicity);
            E.resize(multiplicity);
            E2.resize(multiplicity);
            SI_T.resize(multiplicity);
            T.resize(multiplicity);
            T2.resize(multiplicity);
            MG.resize(multiplicity);
            M.resize(multiplicity);
            Z.resize(multiplicity);
            Indentified.resize(multiplicity);
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
    std::unordered_map<std::string, std::unordered_map<std::string, NPL::EnergyLoss *>> energy_loss;

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
        //Initialization of basic structure
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            fragment->Indentified[ii] = false;
            fragment->Pos[ii] = TVector3((**(data->Mugast)).PosX[ii],
                                         (**(data->Mugast)).PosY[ii],
                                         (**(data->Mugast)).PosZ[ii]);
            fragment->SI_E[ii] = (**(data->Mugast)).DSSD_E[ii];
            fragment->SI_E2[ii] = (**(data->Mugast)).SecondLayer_E[ii];
            fragment->SI_X[ii] = (**(data->Mugast)).DSSD_X[ii];
            fragment->SI_Y[ii] = (**(data->Mugast)).DSSD_Y[ii];
            fragment->SI_T[ii] = (**(data->Mugast)).DSSD_T[ii];
            fragment->T2[ii] = (**(data->Mugast)).SecondLayer_T[ii];
            fragment->MG[ii] = (**(data->Mugast)).TelescopeNumber[ii];
        }
#ifdef VERBOSE_DEBUG
        std::cout << "------------>finished:setting up fragment\n";
#endif

        //Applying (time) re-calibrations
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            if (calibrations_TY[fragment->MG[ii]] == nullptr)
                fragment->T[ii] = fragment->SI_T[ii];
            else
                fragment->T[ii] = calibrations_TY[fragment->MG[ii]]
                                      ->Evaluate(fragment->T[ii], fragment->SI_Y[ii]);
        };
#ifdef VERBOSE_DEBUG
        std::cout << "------------>finished: calibrations\n";
#endif

        //Identification with E TOF
        try {
            for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
                for (const auto &cut_it : particles) {
                    if (with_cuts && 
                            cut_type["E_TOF"].at("E_TOF_" +
                                                    cut_it + 
                                                    "_MG" + 
                                                    std::to_string(fragment->MG[ii]))
                                                        ->IsInside(fragment->SI_E[ii], fragment->T[ii])) {

                        if (fragment->Indentified[ii]) 
                            throw std::runtime_error("Overlapping MUGAST E TOF gates");

                        fragment->M[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("m") + 1,
                                                                  cut_it.find_first_of("_") - cut_it.find_first_of("m") - 1));

                        fragment->Z[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("z") + 1,
                                                                  cut_it.find_first_of("_") - cut_it.find_first_of("z") - 1));

                        fragment->Indentified[ii] = true;
                    }
                }
                if (!fragment->Indentified[ii]) {
                    fragment->M[ii] = 0;
                    fragment->Z[ii] = 0;
                }
            };
        } catch (const std::out_of_range &err) {
            std::cerr << "Mugast cuts not found\n";
            with_cuts = false;
        }

        //Energy reconstruction
        std::unordered_map<std::string, NPL::EnergyLoss *> *ptr_tmp;
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            if (!fragment->Indentified[ii]) {
                fragment->E[ii] = fragment->SI_E[ii];
                continue;
            }
            double tmp_en = 0;

            ptr_tmp = &energy_loss["m" +
                                   std::to_string(fragment->M[ii]) +
                                   "_z" +
                                   std::to_string(fragment->Z[ii])];

            TVector3 hit_direction = fragment->Pos[ii] - target_pos;
            double theta = hit_direction.Angle(TVector3(0, 0, -1));

            //Passivation layer
            tmp_en = (*ptr_tmp)["al_front"]
                         ->EvaluateInitialEnergy(fragment->SI_E[ii],
                                                 0.4E-3,  //Units in mm!
                                                 hit_direction.Angle((**(data->Mugast)).GetTelescopeNormal(ii)));
            tmp_en = (*ptr_tmp)["ice_front"]
                         ->EvaluateInitialEnergy(tmp_en,
                                                 15E-3,  //Units mm!
                                                 fragment->Pos[ii].Angle(TVector3(0, 0, -1)));
            tmp_en = (*ptr_tmp)["havar_front"]
                         ->EvaluateInitialEnergy(tmp_en,
                                                 3.8E-3,
                                                 havar_angle->Evaluate(theta));
            tmp_en = (*ptr_tmp)["he3_front"]
                         ->EvaluateInitialEnergy(tmp_en,
                                                 gas_thickness->Evaluate(theta),
                                                 0.);

            fragment->E[ii] = tmp_en;
        }

        //TODO: compute Ex here
        return true;
    }

    inline int Get_Mult() { return fragment->multiplicity; };
    inline TVector3 *Get_Pos(const int &i) { return &(fragment->Pos[i]); };
    inline int Get_SI_X(const int &i) { return fragment->SI_X[i]; };
    inline int Get_SI_Y(const int &i) { return fragment->SI_Y[i]; };
    inline double Get_SI_E(const int &i) { return fragment->SI_E[i]; };
    inline double Get_T2(const int &i) { return fragment->T2[i]; };
    inline double Get_MG(const int &i) { return fragment->MG[i]; };
    inline double Get_M(const int &i) { return fragment->M[i]; };
    inline double Get_Z(const int &i) { return fragment->Z[i]; };
    inline double Get_E(const int &i) { return fragment->E[i]; };
    inline double Get_T(const int &i) { return fragment->T[i]; }
};

#endif