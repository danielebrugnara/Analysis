#ifndef __REACTIONFRAGMENT_H__
#define __REACTIONFRAGMENT_H__

#include <string>
#include <fstream>

#include <TLorentzVector.h>
#include <TVector3.h>

#include <Units.h>

class ReactionFragment{
public:
    struct FragmentSettings{
        unsigned int    A;   //A
        unsigned int    Z;   //Z
        unsigned int    Q;   //Q
        double          Ek;  //ELab
        double          Ex;  //ELab
        TVector3        Dir; //Dir
        FragmentSettings(   const unsigned int &  A,
                            const unsigned int &  Z,
                            const unsigned int &  Q,
                            const double       &  Ek,
                            const double       &  Ex,
                            const TVector3     &  Dir):
                A(A),
                Z(Z),
                Q(Q),
                Ek(Ek),
                Ex(Ex),
                Dir(Dir){};
    };
    ReactionFragment(FragmentSettings const &);
    ~ReactionFragment()=default;

private:
    ReactionFragment()=default;
    unsigned int A;
    unsigned int Z;
    unsigned int Q;
    double M;
    double Ek;
    double E;
    double Ex;
    double Spin;
    double LifeTime;
    double MassExcess;
    std::string Parity;
    TVector3 Pos;
    std::string Name;
    TLorentzVector P;

public:
    //Setter methods
    inline void SetEx(double const & Ex) {this->Ex = Ex; M+=Ex;};
    inline void SetPos(TVector3 const & Pos) {this->Pos = Pos;};

    //Getter methods
    inline unsigned int Get_A          (){return A;};
    inline unsigned int Get_Z          (){return Z;};
    inline unsigned int Get_Q          (){return Q;};
    inline double       Get_M          (){return M;};
    inline double       Get_Ek         (){return Ek;};
    inline double       Get_E          (){return E;};
    inline double       Get_Ex         (){return Ex;};
    inline double       Get_Spin       (){return Spin;};
    inline double       Get_LifeTime   (){return LifeTime;};
    inline double       Get_MassExcess (){return MassExcess;};
    inline std::string  Get_Parity     (){return Parity;};
    inline TVector3     Get_Pos        (){return Pos;};
    inline std::string  Get_Name       (){return Name;};
    inline TLorentzVector Get_P        (){return P;};

private:
    void GetFromData(const int &, const int &);
    void Extract(std::string);
};

#endif