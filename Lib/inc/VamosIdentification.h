
#ifndef __VAMOSIDENTIFICATION_H__
#define __VAMOSIDENTIFICATION_H__

#include "Identification.h"
#include "Interpolation.h"
#include <TLorentzVector.h>

#include <array>
#include <sstream>

class VamosIdentification : public Identification {
   public:
    VamosIdentification();
    ~VamosIdentification();

    struct Data {
        const double *IC;
        //const double *MW_N;
        const double Path;
        const double Brho;
        const double Xf;
        const unsigned long AGAVA_VAMOSTS;
        const double T_FPMW_CATS2_C;
        Data(   const float * IC, 
                const double Path,
                const double Brho,
                const double Xf,
                const unsigned long AGAVA_VAMOSTS,
                const double T_FPMW_CATS2_C):
                IC(IC), 
                Path(Path),
                Brho(Brho),
                Xf(Xf),
                AGAVA_VAMOSTS(AGAVA_VAMOSTS),
                T_FPMW_CATS2_C(T_FPMW_CATS2_C){};
    };
    Data const* data; 


   private:
    void ReadFPTimeShifts();
    inline double GetShift();
    inline double GetFPTime();
    std::vector<std::pair<double, double>> TimeShifts; //Xf-max, shift
    struct Fragment {
        double En;            //Energy IC
        double D_En;          //DE 1+2
        double D_En2;         //DE 1
        double Path;          //Reconstructed path
        double T;             //Time
        double V;             //Velocity
        double Beta;          //Beta
        double Gamma;         //Gamma
        double M_Q;           //Mass over charge
        double M;             //Mass
        double Charge;        //Charge state
        bool Identified;      //Positive identification
        std::string Nucleus;  //Name in format 47_K
        TLorentzVector p4;    //4 momentum of recoil
        Fragment(): En(0),
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
                    Nucleus(""),
                    p4(){};
    };
    Fragment *fragment;

    Interpolation *FP_time_interpolation;

    const double IC_threashold {0.1};

    public:
    inline void SetData(const Data &);
    inline void Identify();

    //Various computed fragment properties
    inline double Get_En()      {return fragment->En;};    
    inline double Get_D_En()    {return fragment->D_En;};  
    inline double Get_D_En2()   {return fragment->D_En2;}; 
    inline double Get_Path()    {return fragment->Path;};  
    inline double Get_T()       {return fragment->T;};     
    inline double Get_V()       {return fragment->V;};     
    inline double Get_Beta()    {return fragment->Beta;};  
    inline double Get_Gamma()   {return fragment->Gamma;}; 
    inline double Get_M_Q()     {return fragment->M_Q;};   
    inline double Get_M()       {return fragment->M;};     
    inline double Get_Charge()  {return fragment->Charge;};

};

#endif
