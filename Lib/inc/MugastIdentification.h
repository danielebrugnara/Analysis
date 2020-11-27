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
#include "MugastData.h"
#include "Must2Data.h"
#include "CatsData.h"

#include "TMugastPhysics.h"
#include "TMust2Physics.h"
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
        TTreeReaderValue<TMust2Physics> *Must2;
        TTreeReaderValue<TCATSPhysics> *Cats;
        TTreeReaderValue<float> *TW;
        int VAMOS_id_M;
        int VAMOS_id_Z;
        Data(TTreeReaderValue<TMugastPhysics> *Mugast,
             TTreeReaderValue<TMust2Physics> *Must2,
             TTreeReaderValue<TCATSPhysics> *Cats,
             TTreeReaderValue<float> *TW,
             int VAMOS_id_M,
             int VAMOS_id_Z) :  Mugast(Mugast),
                                Must2(Must2),
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


private: //Data input and output
    Data const *data;
    MugastData fragment;
    Must2Data fragment_must;
    CatsData fragment_cats;

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
    inline unsigned int Get_Mult()             const { return fragment.multiplicity; };
    inline TVector3*    Get_Pos(const int &i)        { return &(fragment.Pos[i]); };
    inline int          Get_SI_X(const int &i) const { return fragment.SI_X[i]; };
    inline int          Get_SI_Y(const int &i) const { return fragment.SI_Y[i]; };
    inline double       Get_SI_E(const int &i) const { return fragment.SI_E[i]; };
    inline double       Get_T2(const int &i)   const { return fragment.T2[i]; };
    inline double       Get_MG(const int &i)   const { return fragment.MG[i]; };
    inline double       Get_M(const int &i)    const { return fragment.M[i]; };
    inline double       Get_Z(const int &i)    const { return fragment.Z[i]; };
    inline double       Get_E(const int &i)    const { return fragment.E[i]; };
    inline double       Get_Ex(const int &i)   const { return fragment.Ex[i]; };
    inline double       Get_T(const int &i)    const { return fragment.T[i]; }
    inline std::string  Get_Particle(const int &i) const { return fragment.Particle[i]; }
    inline double       Get_ThetaCM(const int &i)  const {return fragment.E_CM[i];};
    inline TVector3*    Get_EmissionDirection(const int &i)
        { return &(fragment.EmissionDirection[i]); };
    inline double       Get_ThetaLab(const int &i) const
        {return fragment.EmissionDirection[i].Angle(TVector3(0, 0, 1));};
    inline double       Get_Phi(const int &i) const
        {return fragment.EmissionDirection[i].Phi();};

    inline MugastData&  Get_MG_Data()                   {return  fragment;}
    inline Must2Data&  Get_MM_Data()                    {return  fragment_must;}
    inline CatsData&  Get_Cats_Data()                   {return  fragment_cats;}

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
        this->data = data;

        //The following constructs at the same address
        fragment.~MugastData();
        new(&fragment) MugastData((**(data->Mugast)).DSSD_E.size());
        fragment_must.~Must2Data();
        new(&fragment_must) Must2Data((**(data->Must2)).EventMultiplicity);
        fragment_cats.~CatsData();
        new(&fragment_cats) CatsData((**(data->Cats)).PositionX.size());
    };
};
