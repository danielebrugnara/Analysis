#pragma once

#include <array>
#include <unordered_map>
#include <map>

#include <TTreeReaderValue.h>
#include <TVector3.h>

#include "Calibration.h"
#include "Identification.h"
#include "Minimizer.h"
#include "Interpolation.h"
#include "EnergyLoss.h"
#include "ReactionReconstruction.h"
#include "ReactionFragment.h"

#include "TMugastPhysics.h"
#include "TCATSPhysics.h"

class MugastIdentification : public Identification
{
public:
    MugastIdentification();
    ~MugastIdentification();
    bool Initialize(const double &, const TVector3 &); //Beam energy in MeV

public: //Data exchange  types
    struct Data
    {
        TTreeReaderValue<TMugastPhysics> *Mugast;
        TTreeReaderValue<TCATSPhysics> *Cats;
        TTreeReaderValue<float> *TW;
        int VAMOS_id_M;
        int VAMOS_id_Z;
        Data(TTreeReaderValue<TMugastPhysics> *Mugast,
             TTreeReaderValue<TCATSPhysics> *Cats,
             TTreeReaderValue<float> *TW,
             int VAMOS_id_M,
             int VAMOS_id_Z) : Mugast(Mugast),
                               Cats(Cats),
                               TW(TW),
                               VAMOS_id_M(VAMOS_id_M),
                               VAMOS_id_Z(VAMOS_id_Z){};
    };

private: //Constants used internally
    static constexpr int n_detectors{6};
    static constexpr int n_strips{128};
    const int charge_state_interpolation{16}; //Cannot be made static constexpr because it is passed as reference
    static constexpr double gradient_descent_normalization{1.E3};
    static constexpr double ice_percentage_second {0.85};
    const double average_beam_thickness = 3.723*UNITS::mm; //500 Torr

public:
    std::array<std::string, 3> light_particles;
    std::array<std::string, 2> strips;
    std::array<int, n_detectors> cuts_MG;

private: //Variables used internally
    //std::array<int, n_detectors> cuts_MG;
    std::array<int, 3> cuts_M;
    std::array<int, 2> cuts_Z;
    //std::array<std::string, 3> light_particles;
    std::array<std::string, 5> fragments;
    //std::array<std::string, 2> strips;

    double beam_energy{};
    double initial_beam_energy{};
    double final_beam_energy{};
    double beam_energy_match_threashold{};
    double brho{};
    TVector3 target_pos;

    std::unordered_map<std::string, unique_ptr<ReactionReconstruction2body>> reaction;
    std::unordered_map<std::string, unique_ptr<ReactionReconstruction2body>>::iterator reaction_it;
    ReactionFragment* beam_ref{};
    std::vector<std::string> layers;

    //Interpolations
    Interpolation *gas_thickness{};
    Interpolation *gas_thickness_cartesian{};
    Interpolation *havar_angle{};
    Interpolation *TW_Brho_M46_Z18{};
    Interpolation *ice_thickness{};
    Interpolation *TW_vs_ice_thickness{};
    std::vector<std::pair<double, double>> TW_vs_ice;

    //Minimizer for ice thickness estimation
    std::unique_ptr<Minimizer> ice_thickness_minimizer;
    std::unique_ptr<Minimizer> distance_minimizer;

    //Calibrations
    std::unordered_map<int, Calibration *> calibrations_TY;

    //Energy Loss
    std::unordered_map<std::string, std::unordered_map<std::string, EnergyLoss *>> energy_loss;
    std::pair<double, double> current_ice_thickness;//first layer, second layer
    bool use_constant_thickness;

    bool with_cuts;
    double havar_thickness;


private: //Analyzed structure
    struct Fragment
    {
        const unsigned int multiplicity;            //Number of particles detected
        std::vector<TVector3> Pos;                  //3D positions
        std::vector<TVector3> EmissionDirection;    //Depends on the target position
        std::vector<TVector3> TelescopeNormal;      //
        std::vector<int> SI_X;                      //Intrinsic X
        std::vector<int> SI_Y;                      //Intrinsic Y
        std::vector<double> SI_E;                   //Energy deposition measured
        std::vector<double> SI_E2;                  //Second layer energy deposition
        std::vector<double> E;                      //After E-Loss corrections
        std::vector<double> E_CM;                   //After E-Loss corrections
        std::vector<double> Ex;                     //Ex computed with reaction
        std::vector<double> E2;                     //Second layer 
        std::vector<double> SI_T;                   //Un calibrated time
        std::vector<double> T;                      //Calibrated time
        std::vector<double> T2;                     //
        std::vector<double> MG;                     //
        std::vector<int> M;                         //Mass number
        std::vector<int> Z;                         //Z number
        std::vector<bool> Indentified;              //
        std::vector<std::string> Particle;          //
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
            E_CM.resize(multiplicity);
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

private: //Data input and output
    Data const *data;
    Fragment *fragment;

private: //Internal initialization methods
    bool InitializeCuts();
    bool InitializeCalibration();
    bool InitializeELoss();
    bool InitializeInterpolations();

private: //Functions used internally during analysis
    void    IdentifyIceThickness();
    double  ComputeDistanceInGas(const TVector3&, const TVector3&);
    double  InitialBeamEnergy(double, double);
    double  MiddleTargetBeamEnergy(double);

public: //Functions called by selector
    bool Identify();
    //Getter methods
    inline unsigned int Get_Mult()              { return fragment->multiplicity; };
    inline TVector3*    Get_Pos(const int &i)   { return &(fragment->Pos[i]); };
    inline int          Get_SI_X(const int &i)  { return fragment->SI_X[i]; };
    inline int          Get_SI_Y(const int &i)  { return fragment->SI_Y[i]; };
    inline double       Get_SI_E(const int &i)  { return fragment->SI_E[i]; };
    inline double       Get_T2(const int &i)    { return fragment->T2[i]; };
    inline double       Get_MG(const int &i)    { return fragment->MG[i]; };
    inline double       Get_M(const int &i)     { return fragment->M[i]; };
    inline double       Get_Z(const int &i)     { return fragment->Z[i]; };
    inline double       Get_E(const int &i)     { return fragment->E[i]; };
    inline double       Get_Ex(const int &i)    { return fragment->Ex[i]; };
    inline double       Get_T(const int &i)     { return fragment->T[i]; }
    inline std::string  Get_Particle(const int &i)  { return fragment->Particle[i]; }
    inline double       Get_ThetaCM(const int &i)   {return fragment->E_CM[i];};
    inline TVector3*    Get_EmissionDirection(const int &i)
        { return &(fragment->EmissionDirection[i]); };
    inline double       Get_ThetaLab(const int &i)
        {return fragment->EmissionDirection[i].Angle(TVector3(0, 0, 1));};

    //Other methods
    std::vector<std::pair<double, double>> GetTWvsIce();

    inline void StoreTWvsIce(){
        if (TW_vs_ice_thickness == nullptr){
            TW_vs_ice.emplace_back(**(data->TW),
                                   current_ice_thickness.first);
        }
    };

    inline void SetData(Data const *data){
        delete this->data;
        delete this->fragment;
        this->data = data;
        fragment = new Fragment((**(data->Mugast)).DSSD_E.size());
    };
};
