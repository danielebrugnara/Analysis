
#ifndef __MUGASTIDENTIFICATION_H__
#define __MUGASTIDENTIFICATION_H__

#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <array>

#include "Calibration.h"
#include "Identification.h"
#include "Minimizer.h"
#include "Interpolation.h"
#include "EnergyLoss.h"

//TODO: generate internal library for the following headers
#include "NPReaction.h"
#include "TMugastPhysics.h"

class MugastIdentification : public Identification
{
public:
    MugastIdentification();
    ~MugastIdentification();
    bool Initialize(const double &, const TVector3 &); //Beam energy in MeV

    static constexpr int n_detectors{6};
    static constexpr int n_strips{128};

    static constexpr double AMU_TO_MEV{931.49436};

    static constexpr int charge_state_interpolation{16};
    static constexpr double gradient_descent_normalization{1.E3};

    static constexpr double ice_percentage_second = 1.;
    struct Data
    {
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
    std::array<std::string, 3> light_particles;
    std::array<std::string, 5> fragments;
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

    struct Fragment
    {
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
        Fragment(const unsigned int multiplicity)
            : multiplicity(multiplicity)
        {
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
    Interpolation *gas_thickness;
    Interpolation *havar_angle;
    Interpolation *TW_Brho_M46_Z18;
    Interpolation *ice_thickness;
    Interpolation *TW_vs_ice_thickness;
    std::vector<std::pair<double, double>> TW_vs_ice;

    //Minimizer for ice thickness estimation
    Minimizer *ice_thickness_minimizer;

    //Calibrations
    std::unordered_map<int, Calibration *> calibrations_TY;

    //Physics
    const double _TO_MEV{931.4936148};
    std::unordered_map<int, std::unordered_map<int, double>> mass;

    //Energy Loss
    std::unordered_map<std::string, std::unordered_map<std::string, EnergyLoss *>> energy_loss;
    std::pair<double, double> current_ice_thickness;
    bool use_constant_thickness;

    bool with_cuts;
    double havar_thickness;
    bool InitializeCuts();
    bool InitializeCalibration();
    bool InitializeELoss();

public:
    std::vector<std::pair<double, double>> GetTWvsIce();
    inline void StoreTWvsIce()
    {
        if (TW_vs_ice_thickness == nullptr)
        {
            TW_vs_ice.push_back(std::pair<double, double>(**(data->TW),
                                                          current_ice_thickness.first));
        }
    };

    inline void SetData(Data const *data)
    {
        if (this->data != nullptr)
            delete this->data;
        if (this->fragment != nullptr)
            delete this->fragment;
        this->data = data;
        fragment = new Fragment((**(data->Mugast)).DSSD_E.size());
    };

    bool Identify();

private:
    void IdentifyIceThickness();
    static constexpr double average_beam_thickness = 3.776; //500 Torr
    double InitialBeamEnergy(double);
    double MiddleTargetBeamEnergy(double);

public: //Functions called by selector
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
    inline double Get_ThetaLab(const int &i)
    {
        return fragment->EmissionDirection[i].Angle(TVector3(0, 0, 1));
    }
};

#endif
