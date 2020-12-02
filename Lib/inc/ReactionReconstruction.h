#pragma once

//Temporary
#include <iostream>
////
#include <string>

#include <TLorentzVector.h>
#include <TVector3.h>

#include <ReactionFragment.h>

template<typename T>
class ReactionReconstruction{
protected:
    ReactionReconstruction( ReactionFragment::FragmentSettings const &,
                            ReactionFragment::FragmentSettings const &);

    virtual void UpdateVectors();
    ReactionFragment p1;
    ReactionFragment p2;
    TLorentzVector P_TOT;
    bool changed_initial_conditions;
    virtual ~ReactionReconstruction()=default;
    const double precision{1.E-8};
public:
    ReactionReconstruction()=delete;
    virtual void SetBeamEnergy(T const &);
};

//p1 + p2 -> p3 + p4 ## where p2 is at rest in lab
template<typename T>
class ReactionReconstruction2body : protected ReactionReconstruction<T>{
private:
    void CheckVectors();
    virtual void UpdateVectors(const bool&);

public:
    ReactionReconstruction2body()=delete;
    typedef std::array<ReactionFragment::FragmentSettings, 4> ReactionInput2body;
    explicit ReactionReconstruction2body(ReactionInput2body const &);
    virtual ~ReactionReconstruction2body()=default;

    virtual std::string Get_Name() const;
    void ChooseFixed(int const &);
    void ChooseExFixed(int const &);
    ReactionFragment& GetFreeFragment();
    ReactionFragment& GetFixedFragment();
    ReactionFragment& GetExFreeFragment();
    ReactionFragment& GetExFixedFragment();

    T Get_ThetaMax();

    void SetBeamEnergy(T const & E) override;

    void Set_E              (T const &);
    void Set_Ek             (T const &);
    T Set_E_cm              (T const &);
    T Set_Ek_cm             (T const &);
    void Set_Theta          (T const &, bool const &);
    void Set_Theta_cm       (T const &);
    void Set_Theta_Phi      (T const &, T const &, bool const &);
    void Set_Theta_Phi_cm   (T const &, T const &);
    T Set_E_Theta           (T const &, T const &, const bool& = true);
    T Set_Ek_Theta          (T const &, T const &, const bool& = true);
    T Set_E_Theta_cm        (T const &, T const &);
    T Set_Ek_Theta_cm       (T const &, T const &);
    T Set_E_Theta_Phi       (T const &, T const &, T const &);
    T Set_Ek_Theta_Phi      (T const &, T const &, T const &);
    T Set_E_Theta_Phi_cm    (T const &, T const &, T const &);
    T Set_Ek_Theta_Phi_cm   (T const &, T const &, T const &);
    void Set_Beta           (TVector3 const &);
    void Set_P              (T const &);

    ReactionFragment & GetReactionFragment(const int & nr) {
        switch (nr) {
            case 1: return this->p1;
            case 2: return this->p2;
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
template<typename T>
class ReactionReconstruction3body : protected ReactionReconstruction2body<T>{
private:

public:
    ReactionReconstruction3body()=delete;
    ReactionReconstruction3body(std::array<ReactionFragment::FragmentSettings, 4> const &,
                                ReactionFragment::FragmentSettings const &,
                                ReactionFragment::FragmentSettings const &);
    virtual ~ReactionReconstruction3body()=default;

    std::string Get_Name() const override;
protected:
    ReactionFragment p5;
    ReactionFragment p6;
    bool CheckConsistency() const;
};