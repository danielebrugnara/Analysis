
#ifndef __MUGASTIDENTIFICATION_H__
#define __MUGASTIDENTIFICATION_H__

#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <array>

#include "Calibration.h"
#include "Identification.h"
#include "Minimizer.h"
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

    static constexpr int n_detectors	{6};
    static constexpr int n_strips 	{128};

    static constexpr double AMU_TO_MEV	{931.49436};

    static constexpr int charge_state_interpolation{16};

    struct Data {
        TTreeReaderValue<TMugastPhysics> *Mugast;
        TTreeReaderValue<float> *TW;
        int VAMOS_id_M;
        int VAMOS_id_Z;
        Data(TTreeReaderValue<TMugastPhysics> *Mugast,
             TTreeReaderValue<float> *TW,
             int VAMOS_id_M,
             int VAMOS_id_Z) : Mugast(Mugast),
                                TW(TW),
                                VAMOS_id_M(VAMOS_id_M),
                                VAMOS_id_Z(VAMOS_id_Z){};
    };

    std::array<int, n_detectors> cuts_MG;
    std::array<int, 3> cuts_M;
    std::array<int, 2> cuts_Z;
    std::array<std::string, 3> particles;
    std::array<std::string, 2> strips;

   private:
    double beam_energy;
    double initial_beam_energy;
    double final_beam_energy;
    double beam_energy_match_threashold;
    double brho;
    TVector3 target_pos;
    std::unordered_map<std::string, NPL::Reaction *> reaction;
    std::unordered_map<std::string, NPL::Reaction *>::iterator reaction_it;
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
        std::vector<double> Ex;
        std::vector<double> E2;
        std::vector<double> SI_T;
        std::vector<double> T;
        std::vector<double> T2;
        std::vector<double> MG;
        std::vector<int> M;
        std::vector<int> Z;
        std::vector<bool> Indentified;
        std::vector<std::string> Particle;
        Fragment(const unsigned int multiplicity) : multiplicity(multiplicity) {
            Pos.resize(multiplicity);
            EmissionDirection.resize(multiplicity);
            TelescopeNormal.resize(multiplicity);
            SI_X.resize(multiplicity);
            SI_Y.resize(multiplicity);
            SI_E.resize(multiplicity);
            SI_E2.resize(multiplicity);
            E.resize(multiplicity);
            Ex.resize(multiplicity);
            E2.resize(multiplicity);
            SI_T.resize(multiplicity);
            T.resize(multiplicity);
            T2.resize(multiplicity);
            MG.resize(multiplicity);
            M.resize(multiplicity);
            Z.resize(multiplicity);
            Indentified.resize(multiplicity);
            Particle.resize(multiplicity);
        };
    };
    Data const *data;
    Fragment *fragment;

    //Interpolations
    Interpolation* gas_thickness;
    Interpolation* havar_angle;
    Interpolation* TW_Brho_M46_Z18;

    //Minimizer for ice thickness estimation
    Minimizer* ice_thickness_minimizer;

    //Calibrations
    std::unordered_map<int, Calibration *> calibrations_TY;

    //Physics
    const double _TO_MEV{931.4936148};
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

        //Evaluate Ice thickness
        IdentifyIceThickness();

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
            fragment->Particle[ii] = "NONE";
            for (const auto &cut_it : particles) {
                if (with_cuts) {
                    try {
                        tmp_cut = cut_type["E_TOF"].at("E_TOF_" + cut_it + "_MG" +
                                                       std::to_string(static_cast<int>(fragment->MG[ii])));
                    } catch (const std::out_of_range &err) {
                        std::cerr << "Mugast cuts not found : "
                                  << "E_TOF_" + cut_it + "_MG"
                                  << std::to_string(static_cast<int>(fragment->MG[ii]))
                                  << std::endl;
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
                        fragment->Particle[ii] = cut_it;
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

        //Energy reconstruction
        std::unordered_map<std::string, NPL::EnergyLoss *> *ptr_tmp;
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            //if (!fragment->Indentified[ii]) {
            fragment->EmissionDirection[ii] = fragment->Pos[ii] - target_pos;
            if (!fragment->Indentified[ii] || fragment->M[ii] == 4) {  //TODO: fix to include alphas
                fragment->E[ii] = fragment->SI_E[ii];
                fragment->Ex[ii] = 0;
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
                                                 15E-3,  //Units mm!
                                                 //current_ice_thickness,
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

            if ((reaction_it = reaction.find("M" + std::to_string(data->VAMOS_id_M) +
                                            "_Z" + std::to_string(data->VAMOS_id_Z) +
                                            fragment->Particle[ii])) != reaction.end()) {
                fragment->Ex[ii] = reaction_it->second
                                                ->ReconstructRelativistic(fragment->E[ii],
                                                                            Get_ThetaLab(ii));

            } else {
                fragment->Ex[ii] = 0;
            }
        }

        //Ex

#ifdef VERBOSE_DEBUG
        std::cout << "------------>finished : eloss calculations\n";
#endif

        //TODO: compute Ex here
        return true;
    }
	private:
    inline void IdentifyIceThickness(){
        brho = TW_Brho_M46_Z18->Evaluate(**(data->TW));
        final_beam_energy = sqrt(pow(brho/3.3356E-3*charge_state_interpolation,2)+pow(mass[46][18], 2))-(mass[46][18]);
	    initial_beam_energy = InitialBeamEnergy(final_beam_energy);

    //    current_ice_thickness = energy_loss["beam"]["ice_back"]->EvaluateMaterialThickness(final_beam_energy,
    //                                                                initial_beam_energy,
    //                                                                30,
    //                                                                1E-4);
    //    current_ice_thickness = current_ice_thickness/2; 
    //    std::cout << "current ice thickness : " << current_ice_thickness << std::endl;
    //    std::cout << "initial energy : " << initial_beam_energy << std::endl;
    //    std::cout << "brho : " << brho << std::endl;
    //    std::cout << "final energy : " << final_beam_energy << std::endl;
        //delete ice_thickness_minimizer; //TODO : understand why this causes segfault
        if (abs(beam_energy-initial_beam_energy)>beam_energy_match_threashold){
            ice_thickness_minimizer = new Minimizer((beam_energy-initial_beam_energy)/1.E3,
                                                    current_ice_thickness,
                                                    7E-7,
                                                    0.002,
                                                    beam_energy_match_threashold/1.E3 , 
                                                    200, 
                                                    1);
            for (int ii=0; ii<100;++ii){
                current_ice_thickness = 
                    ice_thickness_minimizer
                        ->PerformStep((beam_energy - InitialBeamEnergy(final_beam_energy))/1.E3);
            }
            for (const auto & it: reaction){
                it.second->SetBeamEnergy(MiddleTargetBeamEnergy(final_beam_energy));
//                std::cout << "Beam Energy : " << MiddleTargetBeamEnergy(final_beam_energy)<< std::endl;
//                std::cout << "Ice thickness : " << current_ice_thickness << std::endl;
            }
            //std::cout << "Ice thickness : " <<current_ice_thickness << std::endl;
        }
    }

	inline double InitialBeamEnergy(double beam_energy_from_brho){
//        std::cout << "final energy: " << beam_energy_from_brho << std::endl;
        beam_energy_from_brho =
            energy_loss["beam"]["ice_back"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        current_ice_thickness,
                                        0);
//        std::cout << "before ice: " << beam_energy_from_brho << std::endl;

        beam_energy_from_brho =
            energy_loss["beam"]["havar_back"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        havar_thickness,
                                        0);
//        std::cout << "before havar: " << beam_energy_from_brho << std::endl;

        beam_energy_from_brho =
            energy_loss["beam"]["he3_front"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        2 * gas_thickness->Evaluate(0.),
                                        0);
//        std::cout << "before he("<<2* gas_thickness->Evaluate(0.)<<") : " << beam_energy_from_brho << std::endl;

        beam_energy_from_brho =
            energy_loss["beam"]["havar_front"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        havar_thickness,
                                        0);
//        std::cout << "before havar: " << beam_energy_from_brho << std::endl;

        beam_energy_from_brho =
            energy_loss["beam"]["ice_front"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        current_ice_thickness,
                                        0);
//        std::cout << "before ice ("<<current_ice_thickness<<"): " << beam_energy_from_brho << std::endl;
//        std::cout << "-------->Actual beam energy should be: "<< beam_energy <<std::endl;
	    return beam_energy_from_brho;
	}

	inline double MiddleTargetBeamEnergy(double beam_energy_from_brho){
        beam_energy_from_brho =
            energy_loss["beam"]["ice_back"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        current_ice_thickness,
                                        0);

        beam_energy_from_brho =
            energy_loss["beam"]["havar_back"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        havar_thickness,
                                        0);

        beam_energy_from_brho =
            energy_loss["beam"]["he3_front"]
                ->EvaluateInitialEnergy(beam_energy_from_brho,
                                        gas_thickness->Evaluate(0.),
                                        0);
	    return beam_energy_from_brho;
	}

	public:

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
    inline double Get_Ex(const int &i) { return fragment->Ex[i]; };
    inline double Get_T(const int &i) { return fragment->T[i]; }
    inline std::string Get_Particle(const int &i) { return fragment->Particle[i]; }
    inline double Get_ThetaLab(const int &i) {
        return fragment->EmissionDirection[i].Angle(TVector3(0, 0, 1));
    }
};

#endif
