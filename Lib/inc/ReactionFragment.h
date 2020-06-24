#ifndef __REACTIONFRAGMENT_H__
#define __REACTIONFRAGMENT_H__

#include <string>
#include <fstream>

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
        double          Ek;  //ELab
        double          Ex;  //ELab
        FragmentSettings(   const unsigned int &  A,
                            const unsigned int &  Z,
                            const unsigned int &  Q,
                            const double       &  Ek,
                            const double       &  Ex):
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
    double Spin;
    double LifeTime;
    double MassExcess;
    std::string Parity;
    std::string Name;

    //Can be changed
    double M;
    double M2;
    double Ex;
    double Invariant;
    struct Properties {
        double Ek;  //primary
        double E;
        TVector3 Pos;
        //TVector3 Beta;
        TLorentzVector P;
        bool modified;
        std::pair <bool, TLorentzVector> P_tot; //Primary if present
        Properties( const double & Ek):
            Ek(Ek),
            E(0),
            modified(false){}
    };
    Properties Lab, Cm;
    Properties & GetModified();
    Properties & GetNotModified();
    bool Fixed;
    bool ExFixed;
    bool preserve_pcm;

public:
    //Setter methods
    void Set_Ex             (double const &);
    void Set_Ek             (double const &);
    void Set_Ek_cm          (double const &);
    void Set_E              (double const &);
    void Set_E_cm           (double const &);
    void Set_Pos            (TVector3 const &);
    void Set_Pos_cm         (TVector3 const &);
    void Set_Dir            (TVector3 const &);
    void Set_Dir_cm         (TVector3 const &);
    void Set_Beta           (TVector3 const &);
    void Set_Beta_cm        (TVector3 const &);
    void Set_P              (TLorentzVector const &);
    void Set_P_cm           (TLorentzVector const &);
    void Set_P              (TVector3 const &);
    void Set_P_cm           (TVector3 const &);
    void Set_P              (double const &);
    void Set_P_cm           (double const &);
    void Set_Theta          (double const &);
    void Set_Theta_cm       (double const &);
    void Set_Theta_Phi      (double const &, double const &);
    void Set_Theta_Phi_cm   (double const &, double const &);
    void Set_E_Theta        (double const &, double const &);
    void Set_E_Theta_cm     (double const &, double const &);
    void Set_E_Theta_Phi    (double const &, double const &, double const &);
    void Set_E_Theta_Phi_cm (double const &, double const &, double const &);
    void Set_P_tot          (TLorentzVector const &);
    void Set_P_tot_cm       (TLorentzVector const &);
    void Set_Fixed          (bool const &);
    void Set_ExFixed        (bool const &);
    inline void Set_Invariant (double const & Inv){Invariant=Inv;};
    bool Check_Consistent   ();

    //Getter methods
    inline unsigned int Get_A          () const {return A;}            ;
    inline unsigned int Get_Z          () const {return Z;}            ;
    inline unsigned int Get_Q          () const {return Q;}            ;
    inline double       Get_M          () const {return M;}            ;
    inline double       Get_M_GS       () const {return M-Ex;}         ;
    inline double       Get_M2         () const {return M2;}           ;
    inline double       Get_Ek         () const {return Lab.Ek;}       ;
    inline double       Get_Ek_cm      () const {return Cm.Ek;}        ;
    inline double       Get_E          () const {return Lab.E;}        ;
    inline double       Get_E_cm       () const {return Cm.E;}         ;
    inline double       Get_Ex         () const {return Ex;}           ;
    inline double       Get_Invariant  () const {return Invariant;}    ;
    inline double       Get_Spin       () const {return Spin;}         ;
    inline double       Get_LifeTime   () const {return LifeTime;}     ;
    inline double       Get_MassExcess () const {return MassExcess;}   ;
    inline std::string  Get_Parity     () const {return Parity;}       ;
    inline TVector3     Get_Pos        () const {return Lab.Pos;}      ;
    inline TVector3     Get_Pos_cm     () const {return Cm.Pos;}       ;
    inline std::string  Get_Name       () const {return Name;}         ;
    inline TLorentzVector Get_P        () const {return Lab.P;}        ;
    inline TLorentzVector Get_P_cm     () const {return Cm.P;}         ;
    inline bool         Is_Fixed       () const {return Fixed;}        ;
    inline bool         Is_ExFixed     () const {return ExFixed;}      ;

    double Compute_Ex       (TLorentzVector const &);
private:
    void GetFromData(const int &, const int &);
    void Extract(std::string const &);
    void ComputeChanges();
    const double precision;
};

#endif