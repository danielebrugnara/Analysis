#pragma once

#include <string>
#include <fstream>
#include <cctype>

#include <TLorentzVector.h>

#include <TVector3.h>

#include <Units.h>

class ReactionFragment{
public:
    ReactionFragment()=delete;
    struct FragmentSettings{
        unsigned int    A;   //A
        unsigned int    Z;   //Z
        unsigned int    Q;   //Q
        long double          Ek;  //ELab
        long double          Ex;  //ELab
        FragmentSettings(   const unsigned int &  A,
                            const unsigned int &  Z,
                            const unsigned int &  Q,
                            const long double       &  Ek,
                            const long double       &  Ex):
                A(A),
                Z(Z),
                Q(Q),
                Ek(Ek),
                Ex(Ex){};
    };
    explicit ReactionFragment(FragmentSettings const &);
    ~ReactionFragment()=default;

private:
    //Cannot be changed
    unsigned int A;
    unsigned int Z;
    unsigned int Q;
    long double Spin;
    long double LifeTime;
    long double MassExcess;
    std::string Parity;
    std::string Name;

    //Can be changed
    long double M;
    long double M2;
    long double Ex;
    struct Properties {
        TVector3 Pos;
        TLorentzVector P;
        Properties(){}
    };
    std::pair <bool, TVector3> betacm;
    Properties Lab, Cm;
    bool Fixed;
    bool ExFixed;

public:
    //Setter methods
    void Set_Ex             (long double const &)            ; //Keeps p constant
    void Set_Ek             (long double const &)            ;
    void Set_Ek_cm          (long double const &)            ;
    void Set_E              (long double const &)            ;
    void Set_E_cm           (long double const &)            ;
    void Set_Pos            (TVector3 const &)          ;
    void Set_Pos_cm         (TVector3 const &)          ;
    void Set_Dir            (TVector3 const &)          ;
    void Set_Dir_cm         (TVector3 const &)          ;
    void Set_Beta           (TVector3 const &)          ;
    void Set_Beta_cm        (TVector3 const &)          ;
    void Set_P              (TLorentzVector const &)    ;
    void Set_P_cm           (TLorentzVector const &)    ;
    void Set_P              (TVector3 const &)          ;
    void Set_P_cm           (TVector3 const &)          ;
    void Set_P              (long double const &)            ;
    void Set_P_cm           (long double const &)            ;
    void Set_Theta          (long double const &)            ;
    void Set_Theta_cm       (long double const &)            ;
    void Set_Theta_Phi      (long double const &, long double const &);
    void Set_Theta_Phi_cm   (long double const &, long double const &);
    void Set_E_Theta        (long double const &, long double const &);
    void Set_E_Theta_cm     (long double const &, long double const &);
    void Set_E_Theta_Phi    (long double const &, long double const &, long double const &);
    void Set_E_Theta_Phi_cm (long double const &, long double const &, long double const &);
    void Set_betacm         (TVector3 const &)          ;
    void Set_Fixed          (bool const &)              ;
    void Set_ExFixed        (bool const &)              ;
    bool Check_Consistent   ()                          ;

    //Getter methods
    inline unsigned int      Get_A           () const {return A;}            ;
    inline unsigned int      Get_Z           () const {return Z;}            ;
    inline unsigned int      Get_Q           () const {return Q;}            ;
    inline long double       Get_M           () const {return M;}            ;
    inline long double       Get_M_GS        () const {return M-Ex;}         ;
    inline long double       Get_M2          () const {return M2;}           ;
    inline long double       Get_Ek          () const {return Lab.P.E()-M;}  ;
    inline long double       Get_Ek_cm       () const {return Cm.P.E()-M;}   ;
    inline long double       Get_E           () const {return Lab.P.E();}    ;
    inline long double       Get_E_cm        () const {return Cm.P.E();}     ;
    inline long double       Get_Ex          () const {return Ex;}           ;
    inline long double       Get_Spin        () const {return Spin;}         ;
    inline long double       Get_LifeTime    () const {return LifeTime;}     ;
    inline long double       Get_MassExcess  () const {return MassExcess;}   ;
    inline std::string       Get_Parity      () const {return Parity;}       ;
    inline TVector3          Get_Pos         () const {return Lab.Pos;}      ;
    inline TVector3          Get_Pos_cm      () const {return Cm.Pos;}       ;
    inline std::string       Get_Name        () const {return Name;}         ;
    inline TLorentzVector    Get_P           () const {return Lab.P;}        ;
    inline TLorentzVector    Get_P_cm        () const {return Cm.P;}         ;
    inline bool              Is_Fixed        () const {return Fixed;}        ;
    inline bool              Is_ExFixed      () const {return ExFixed;}      ;
    inline TVector3          Get_betacm      () const {
        if(betacm.first) return betacm.second;
        else throw std::runtime_error("Betacm was not set!\n");
    }

    inline double             GetBrhoFromEk_Q(const long double& Ek, const int& Q) const{
        return (double)sqrtl(Ek*Ek+2.*Ek*M)*UNITS::e_SI /(UNITS::PHYSICS::c_light/1E6) / (Q * UNITS::e_SI);
    }
    inline double             GetEkFromBrho_Q(const long double& Brho, const int& Q) const{
        return (double)(sqrtl(pow(Brho * UNITS::PHYSICS::c_light/1E6*Q, 2)+M2)-M);
    }

    long double Get_BindingEnergy() const;
    void BoostToLab()                           ;
    void BoostToCm()                            ;
private:
    void GetFromData(const int &, const int &)  ;
    void Extract(std::string const &)           ;
    void UpdateMass()                           ;
    const long double precision;

public:
    friend inline bool operator==(const ReactionFragment& lhs, const ReactionFragment& rhs){
        return &lhs==&rhs;
    }
    friend inline bool operator!=(const ReactionFragment& lhs, const ReactionFragment& rhs){
        return &lhs!=&rhs;
    }
};