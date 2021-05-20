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
    bool initialize(const double &beamEnergy, const TVector3 &targetPos); //Beam energy in MeV

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
             TTreeReaderValue<Float_t> *TW,
             int VAMOS_id_M,
             int VAMOS_id_Z) :  Mugast(Mugast),
                                Must2(Must2),
                                Cats(Cats),
                                TW(TW),
                                VAMOS_id_M(VAMOS_id_M),
                                VAMOS_id_Z(VAMOS_id_Z){};
    };

private: //Constants used internally
    const double targetRadius{8*UNITS::mm};
    static constexpr int nDetectors{6};
    static constexpr int nStrips{128};
    const int chargeStateInterpolation{16}; //Cannot be made static constexpr because it is passed as reference
    //static constexpr double gradient_descent_normalization{1.E3};
    //static constexpr double icePercentageSecond {0.85};
    static constexpr double icePercentageSecond {1.};
    const double averageBeamThickness = 3.7210 * UNITS::mm;    //700 mbar
    //const double averageBeamThickness = 6.39 * UNITS::mm;    //70000000 mbar
    const std::array<std::pair<double, double>, 32> beamPositions = {
            //std::make_pair( 0.5*UNITS::mm,  0.5*UNITS::mm),std::make_pair( 0.5*UNITS::mm, -0.5*UNITS::mm),
            //std::make_pair(-0.5*UNITS::mm, -0.5*UNITS::mm),std::make_pair(-0.5*UNITS::mm,  0.5*UNITS::mm),
            //std::make_pair( 1.0*UNITS::mm,  1.0*UNITS::mm),std::make_pair( 1.0*UNITS::mm, -1.0*UNITS::mm),
            //std::make_pair(-1.0*UNITS::mm, -1.0*UNITS::mm),std::make_pair(-1.0*UNITS::mm,  1.0*UNITS::mm),
            //std::make_pair( 1.5*UNITS::mm,  1.5*UNITS::mm),std::make_pair( 1.5*UNITS::mm, -1.5*UNITS::mm),
            //std::make_pair(-1.5*UNITS::mm, -1.5*UNITS::mm),std::make_pair(-1.5*UNITS::mm,  1.5*UNITS::mm),
            //std::make_pair( 2.0*UNITS::mm,  2.0*UNITS::mm),std::make_pair( 2.0*UNITS::mm, -2.0*UNITS::mm),
            //std::make_pair(-2.0*UNITS::mm, -2.0*UNITS::mm),std::make_pair(-2.0*UNITS::mm,  2.0*UNITS::mm)
            //Diagonal values
            std::make_pair( 1.0*UNITS::mm,  1.0*UNITS::mm),std::make_pair( 1.0*UNITS::mm, -1.0*UNITS::mm),//0  1
            std::make_pair(-1.0*UNITS::mm, -1.0*UNITS::mm),std::make_pair(-1.0*UNITS::mm,  1.0*UNITS::mm),//2  3
            std::make_pair( 2.0*UNITS::mm,  2.0*UNITS::mm),std::make_pair( 2.0*UNITS::mm, -2.0*UNITS::mm),//4  5
            std::make_pair(-2.0*UNITS::mm, -2.0*UNITS::mm),std::make_pair(-2.0*UNITS::mm,  2.0*UNITS::mm),//6  7
            std::make_pair( 3.0*UNITS::mm,  3.0*UNITS::mm),std::make_pair( 3.0*UNITS::mm, -3.0*UNITS::mm),//8  9
            std::make_pair(-3.0*UNITS::mm, -3.0*UNITS::mm),std::make_pair(-3.0*UNITS::mm,  3.0*UNITS::mm),//10 11
            std::make_pair( 4.0*UNITS::mm,  4.0*UNITS::mm),std::make_pair( 4.0*UNITS::mm, -4.0*UNITS::mm),//12 13
            std::make_pair(-4.0*UNITS::mm, -4.0*UNITS::mm),std::make_pair(-4.0*UNITS::mm,  4.0*UNITS::mm),//14 15
            //Plane align values
            std::make_pair( 1.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm,  1.0*UNITS::mm),//16 17
            std::make_pair(-1.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm, -1.0*UNITS::mm),//18 19
            std::make_pair( 2.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm,  2.0*UNITS::mm),//20 21
            std::make_pair(-2.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm, -2.0*UNITS::mm),//22 23
            std::make_pair( 3.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm,  3.0*UNITS::mm),//24 25
            std::make_pair(-3.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm, -3.0*UNITS::mm),//26 27
            std::make_pair( 4.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm,  4.0*UNITS::mm),//28 29
            std::make_pair(-4.0*UNITS::mm,  0.0*UNITS::mm),std::make_pair( 0.0*UNITS::mm, -4.0*UNITS::mm) //30 31
    };
    const std::array<double, 32> focusScale = {  1.10, 1.03, 0.96, 0.89, 0.82, 0.75, 0.68, 0.61, //0 - 7
                                                 0.54, 0.47, 0.40, 0.33, 0.26, 0.19, 0.12, 0.05, //8 - 15
                                                -0.02,-0.09,-0.16,-0.23,-0.30,-0.37,-0.44,-0.51, //16- 23
                                                -0.58,-0.65,-0.72,-0.79,-0.86,-0.93,-1.00,-1.07};//24- 31

public:
    std::array<std::string, 5> lightParticles;
    std::array<int, nDetectors> cutsMg;
    std::array<std::string, 2> strips;

private: //Variables used internally
    double beamEnergy{};
    double initialBeamEnergy{};
    double finalBeamEnergy{};
    double beamEnergyMatchThreashold{};
    double brho{};
    double currentMidTargetBeamEnergy{0};
    TVector3 targetPos;

    std::unordered_map<std::string, unique_ptr<ReactionReconstruction2body<long double>>> reaction;
    std::unordered_map<std::string, unique_ptr<ReactionReconstruction2body<long double>>>::iterator reactionIt;
    ReactionFragment* beamRef{};
    std::vector<std::string> layers;

    //Interpolations
    Interpolation *gasThickness{};
    Interpolation *gasThicknessCartesian{};
    Interpolation *gasThicknessDerivativeCartesian{};
    Interpolation *havarAngle{};
    Interpolation *TW_Brho_M46_Z18{};
    //std::vector<std::pair<double, double>> TW_vs_ice;

    //Minimizer for ice thickness estimation
    std::unique_ptr<Minimizer> iceThicknessMinimizer;
    std::unique_ptr<Minimizer> distanceMinimizer;

    //Calibrations
    std::unordered_map<int, Calibration *> calibrationsTy;

    //Energy Loss
    std::unordered_map<std::string, std::unordered_map<std::string, EnergyLoss *>> energyLoss;
    std::pair<double, double> currentIceThickness;//first layer, second layer
    bool useConstantThickness;

    bool with_cuts;
    double havarThickness;


private: //Data input and output
    Data const *data;
    MugastData fragment;
    Must2Data fragmentMust;
    CatsData fragmentCats;

private: //Internal initialization methods
    bool initializeCuts();
    bool initializeCalibration();
    bool initializeELoss();
    bool initializeInterpolations();

private: //Functions used internally during analysis
    void    identifyIceThickness();
    TVector3 computeCollision(TVector3, const TVector3&);
    //double  computeDistanceInGas(TVector3, const TVector3&);
    double  computeAngleOfIncidence(const TVector3&, const TVector3&);
    double  InitialBeamEnergy(double, double);
    double  middleTargetBeamEnergy(double);
    bool identifyMugast();
    bool identifyMust2();
    bool reconstructEnergy(MugastData&);

    template<class D>
    bool fillInitialData(MugastData&, const D*);

public: //Functions called by selector
    bool identify();
    //Getter methods
    inline unsigned int getMult()                   const { return fragment.multiplicity;   };
    inline const TVector3*    getPos(const int &i)  const { return &(fragment.Pos[i]);      };
    inline int          getSiX(const int &i)        const { return fragment.SI_X[i];        };
    inline int          getSiY(const int &i)        const { return fragment.SI_Y[i];        };
    inline double       getSiE(const int &i)        const { return fragment.SI_E[i];        };
    inline double       getT2(const int &i)         const { return fragment.T2[i];          };
    inline double       getMg(const int &i)         const { return fragment.MG[i];          };
    inline double       getM(const int &i)          const { return fragment.M[i];           };
    inline double       getZ(const int &i)          const { return fragment.Z[i];           };
    inline double       getE(const int &i)          const { return fragment.E[i];           };
    inline double       getEx(const int &i)         const { return fragment.Ex[i];          };
    inline double       getT(const int &i)          const { return fragment.T[i];           };
    inline std::string  getParticle(const int &i)   const { return fragment.Particle[i];    };
    inline double       getThetaCm(const int &i)    const { return fragment.E_CM[i];        };
    inline double       getEBeforeLayer(const int &i)       const { return fragment.EBeforeLayer[i];    };
    inline double       getThicknessInLayer(const int &i)   const { return fragment.ThickessInLayer[i]; };
    inline TVector3*    getEmissionDirection(const int &i)
        { return &(fragment.EmissionDirection[i]); };
    inline double       getThetaLab(const int &i) const
        {return fragment.EmissionDirection[i].Angle(TVector3(0, 0, 1));};
    inline double       getPhi(const int &i) const
        {return fragment.EmissionDirection[i].Phi();};
    inline double       getIceThicknessFront()      const { return currentIceThickness.first;};
    inline double       getIceThicknessBack()       const { return currentIceThickness.second;};
    inline double       getMidTargetBeamEnergy()    const { return currentMidTargetBeamEnergy;};

    inline MugastData&  getMgData()                   { return  fragment;}
    inline Must2Data&   getMmData()                   { return  fragmentMust;}
    inline CatsData&    getCatsData()                 { return  fragmentCats;}

    //Other methods
    //std::vector<std::pair<double, double>> getTWvsIce();

    //inline void storeTWvsIce(){
    //        TW_vs_ice.emplace_back(**(data->TW),
    //                               currentIceThickness.first);
    //};

    inline void setData(Data const *data){
        delete this->data;
        this->data = data;

        //The following constructs at the same address
        fragment.~MugastData();
        new(&fragment) MugastData((**(data->Mugast)).DSSD_E.size(), focusScale.size());
        fragmentMust.~Must2Data();
        new(&fragmentMust) Must2Data((**(data->Must2)).EventMultiplicity, focusScale.size());
        fragmentCats.~CatsData();
        new(&fragmentCats) CatsData((**(data->Cats)).PositionX.size(), focusScale.size());
    };
};
