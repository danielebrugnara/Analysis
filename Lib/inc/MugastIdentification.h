
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
#include "NPReaction.h"

class MugastIdentification : public Identification {
   public:
    MugastIdentification();
    ~MugastIdentification();
    bool Initialize(const double &, const TVector3 &);  //Beam energy in MeV

    static constexpr int n_detectors = 6;
    static constexpr int n_strips = 128;

    struct Data {
        TTreeReaderValue<TMugastPhysics> *Mugast;
        TTreeReaderValue<float> *TW;
        Data(TTreeReaderValue<TMugastPhysics> *Mugast,
                TTreeReaderValue<float> *TW) : Mugast(Mugast),
                                                TW(TW){};
    };
    std::array<int, n_detectors> cuts_MG;
    std::array<int, 3> cuts_M;
    std::array<int, 2> cuts_Z;
    std::array<std::string, 3> particles;
    std::array<std::string, 2> strips;

   private:
    double beam_energy;
    TVector3 target_pos;
    NPL::Reaction *reaction;
    std::vector<std::string> layers;

    struct Fragment {
        const unsigned int multiplicity;
        std::vector<TVector3> Pos;
        std::vector<TVector3> EmissionDirection;
        std::vector<TVector3> TelescopeNormal;
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
            EmissionDirection.resize(multiplicity);
            TelescopeNormal.resize(multiplicity);
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
    Interpolation *TW_Brho_M46_Z18;

    //Calibrations
    std::unordered_map<int, Calibration *> calibrations_TY;

    //Physics
    const double AMU_TO_MEV{931.4936148};
    std::unordered_map<int, std::unordered_map<int, double>> mass;

    //Energy Loss
    std::unordered_map<std::string, std::unordered_map<std::string, NPL::EnergyLoss *>> energy_loss;
    double current_ice_thickness;

    bool with_cuts;
    double havar_thickness;

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
            fragment->TelescopeNormal[ii] = TVector3((**(data->Mugast)).TelescopeNormalX[ii],
                                                     (**(data->Mugast)).TelescopeNormalY[ii],
                                                     (**(data->Mugast)).TelescopeNormalZ[ii]);
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
        TCutG *tmp_cut;
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            for (const auto &cut_it : particles) {
                if (with_cuts) {
                    try {
                        tmp_cut = cut_type["E_TOF"].at("E_TOF_" + cut_it + "_MG" +
                                                        std::to_string(static_cast<int>(fragment->MG[ii])));
                    } catch (const std::out_of_range &err) {
                        std::cerr << "Mugast cuts not found : "
                                  <<"E_TOF_" + cut_it + "_MG"
                                  <<std::to_string(static_cast<int>(fragment->MG[ii]))
                                  <<std::endl;
                        with_cuts = false;
                        continue;
                    }
                    if (tmp_cut->IsInside(fragment->SI_E[ii], fragment->T[ii])) {
                        if (fragment->Indentified[ii])
                            throw std::runtime_error("Overlapping MUGAST E TOF gates :" +
                                                     cut_it +
                                                     "\tMG" +
                                                     fragment->MG[ii] +
                                                     "\n");

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
        }
#ifdef VERBOSE_DEBUG
        std::cout << "------------>finished: searching cuts\n";
#endif

        //Evaluate Ice thickness
        double brho = TW_Brho_M46_Z18->Evaluate(**(data->TW)); 
        double beam_energy = 0;

        beam_energy = 
            energy_loss["beam"]["ice_back"]
                ->EvaluateInitialEnergy(beam_energy,
                                        current_ice_thickness,
                                        0);

        beam_energy = 
            energy_loss["beam"]["havar_back"]
                ->EvaluateInitialEnergy(beam_energy,
                                        havar_thickness,
                                        0);

        beam_energy = 
            energy_loss["beam"]["he3_front"]
                ->EvaluateInitialEnergy(beam_energy,
                                        2*gas_thickness->Evaluate(0.),
                                        0);

        beam_energy = 
            energy_loss["beam"]["havar_back"]
                ->EvaluateInitialEnergy(beam_energy,
                                        havar_thickness,
                                        0);

        beam_energy = 
            energy_loss["beam"]["ice_front"]
                ->EvaluateInitialEnergy(beam_energy,
                                        current_ice_thickness,
                                        0);

        //Energy reconstruction
        std::unordered_map<std::string, NPL::EnergyLoss *> *ptr_tmp;
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            //if (!fragment->Indentified[ii]) {
            fragment->EmissionDirection[ii] = fragment->Pos[ii] - target_pos;
            if (!fragment->Indentified[ii] || fragment->M[ii] == 4) {  //TODO: fix to include alphas
                fragment->E[ii] = fragment->SI_E[ii];
                continue;
            }
            double tmp_en = 0;

            ptr_tmp = &energy_loss["m" +
                                   std::to_string(fragment->M[ii]) +
                                   "_z" +
                                   std::to_string(fragment->Z[ii])];

            double theta = fragment->EmissionDirection[ii].Angle(TVector3(0, 0, -1));

            //Passivation layer
            tmp_en = (*ptr_tmp)["al_front"]
                         ->EvaluateInitialEnergy(fragment->SI_E[ii],
                                                 0.4E-3,  //Units in mm!
                                                 fragment->EmissionDirection[ii]
                                                            .Angle(fragment->TelescopeNormal[ii]));
            tmp_en = (*ptr_tmp)["ice_front"]
                         ->EvaluateInitialEnergy(tmp_en,
                                                 //15E-3,  //Units mm!
                                                 current_ice_thickness,
                                                 havar_angle->Evaluate(theta));
            tmp_en = (*ptr_tmp)["havar_front"]
                         ->EvaluateInitialEnergy(tmp_en,
                                                 //3.8E-3,
                                                 havar_thickness,
                                                 havar_angle->Evaluate(theta));
            tmp_en = (*ptr_tmp)["he3_front"]
                         ->EvaluateInitialEnergy(tmp_en,
                                                 gas_thickness->Evaluate(theta),
                                                 0.);

            fragment->E[ii] = tmp_en;
        }

        //Ex

#ifdef VERBOSE_DEBUG
        std::cout << "------------>finished : eloss calculations\n";
#endif

        //TODO: compute Ex here
        return true;
    }

    inline int Get_Mult() { return fragment->multiplicity; };
    inline TVector3 *Get_Pos(const int &i) { return &(fragment->Pos[i]); };
    inline TVector3 *Get_EmissionDirection(const int &i) { return &(fragment->EmissionDirection[i]); };
    inline int Get_SI_X(const int &i) { return fragment->SI_X[i]; };
    inline int Get_SI_Y(const int &i) { return fragment->SI_Y[i]; };
    inline double Get_SI_E(const int &i) { return fragment->SI_E[i]; };
    inline double Get_T2(const int &i) { return fragment->T2[i]; };
    inline double Get_MG(const int &i) { return fragment->MG[i]; };
    inline double Get_M(const int &i) { return fragment->M[i]; };
    inline double Get_Z(const int &i) { return fragment->Z[i]; };
    inline double Get_E(const int &i) { return fragment->E[i]; };
    inline double Get_T(const int &i) { return fragment->T[i]; }
    inline double Get_ThetaLab(const int &i) { 
        return fragment->EmissionDirection[i].Angle(TVector3(0, 0, 1)); 
    }
};

#endif
