#ifndef __REACTIONRECONSTRUCTION_H__
#define __REACTIONRECONSTRUCTION_H__

//Temporary
#include <iostream>
////
#include <string>

#include <TLorentzVector.h>
#include <TVector3.h>

#include <ReactionFragment.h>

class ReactionReconstruction{
protected:
    ReactionReconstruction( ReactionFragment::FragmentSettings const &,
                            ReactionFragment::FragmentSettings const &);

    virtual void UpdateVectors();
    ReactionFragment p1;
    ReactionFragment p2;
    TLorentzVector P_TOT;
    bool changed_initial_conditions;
    ~ReactionReconstruction()=default;
    const double precision{1.E-8};
public:
    ReactionReconstruction()=delete;
    virtual void SetBeamEnergy(long double const &);
};

//p1 + p2 -> p3 + p4 ## where p2 is at rest in lab
class ReactionReconstruction2body : protected ReactionReconstruction{
private:
    void CheckVectors();
    virtual void UpdateVectors(const bool&);

public:
    ReactionReconstruction2body()=delete;
    typedef std::array<ReactionFragment::FragmentSettings, 4> ReactionInput2body;
    explicit ReactionReconstruction2body(ReactionInput2body const &);
    ~ReactionReconstruction2body()=default;

    virtual std::string Get_Name() const;
    void ChooseFixed(int const &);
    void ChooseExFixed(int const &);
    ReactionFragment& GetFreeFragment();
    ReactionFragment& GetFixedFragment();
    ReactionFragment& GetExFreeFragment();
    ReactionFragment& GetExFixedFragment();

    long double Get_ThetaMax();

    void SetBeamEnergy(long double const & E) override;

    void Set_E              (long double const &);
    void Set_Ek             (long double const &);
    long double Set_E_cm         (long double const &);
    long double Set_Ek_cm         (long double const &);
    void Set_Theta          (long double const &, bool const &);
    void Set_Theta_cm       (long double const &);
    void Set_Theta_Phi      (long double const &, long double const &, bool const &);
    void Set_Theta_Phi_cm   (long double const &, long double const &);
    long double Set_E_Theta      (long double const &, long double const &, const bool& = true);
    long double Set_Ek_Theta      (long double const &, long double const &, const bool&);
    long double Set_E_Theta_cm   (long double const &, long double const &);
    long double Set_Ek_Theta_cm   (long double const &, long double const &);
    long double Set_E_Theta_Phi      (long double const &, long double const &, long double const &);
    long double Set_Ek_Theta_Phi      (long double const &, long double const &, long double const &);
    long double Set_E_Theta_Phi_cm   (long double const &, long double const &, long double const &);
    long double Set_Ek_Theta_Phi_cm   (long double const &, long double const &, long double const &);
    void Set_Beta           (TVector3 const &);
    void Set_P              (long double const &);

    ReactionFragment & GetReactionFragment(const int & nr) {
        switch (nr) {
            case 1: return p1;
            case 2: return p2;
            case 3: return p3;
            case 4: return p4;
            default: throw std::runtime_error("Reaction fragment chosen not existent\n");
        }
    };
protected:
    ReactionFragment p3;
    ReactionFragment p4;
    bool CheckConsistency() const;
};

//p1 + p2 -> p3 + p4* -> p3 + p5 + p6
class ReactionReconstruction3body : protected ReactionReconstruction2body{
private:

public:
    ReactionReconstruction3body()=delete;
    ReactionReconstruction3body(ReactionInput2body const &,
                                ReactionFragment::FragmentSettings const &,
                                ReactionFragment::FragmentSettings const &);
    ~ReactionReconstruction3body()=default;

    std::string Get_Name() const override;
protected:
    ReactionFragment p5;
    ReactionFragment p6;
    bool CheckConsistency() const;
};

#endif