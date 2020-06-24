#ifndef __REACTIONRECONSTRUCTION_H__
#define __REACTIONRECONSTRUCTION_H__

#include <string>

#include <TLorentzVector.h>
#include <TVector3.h>

#include <ReactionFragment.h>

class ReactionReconstruction{
protected:
    ReactionReconstruction( ReactionFragment::FragmentSettings const &,
                            ReactionFragment::FragmentSettings const &);
    ~ReactionReconstruction()=default;
    void UpdateVectors();
    ReactionFragment p1;
    ReactionFragment p2;
    TLorentzVector P_TOT;
    bool changed_initial_conditions;

private:

public:
    //Nothing public
    ReactionReconstruction()=delete;

    virtual void SetBeamEnergy(double const &);
    void Get_P_CM() const ;
};

//p1 + p2 -> p3 + p4 ## where p2 is at rest
class ReactionReconstruction2body : protected ReactionReconstruction{
private:
    void UpdateVectors();

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

    double Get_ThetaMax();

    inline void SetBeamEnergy(double const & E) override;

    void Set_E              (double const &);
    void Set_Ek             (double const &);
    double Set_E_cm         (double const &);
    void Set_Theta          (double const &, bool const &);
    void Set_Theta_cm       (double const &);
    void Set_Theta_Phi      (double const &, double const &, bool const &);
    void Set_Theta_Phi_cm   (double const &, double const &);
    double Set_E_Theta      (double const &, double const &, bool const &);
    double Set_E_Theta_cm   (double const &, double const &);
    double Set_E_Theta_Phi      (double const &, double const &, double const &);
    double Set_E_Theta_Phi_cm   (double const &, double const &, double const &);
    void Set_Dir            (TVector3 const &);
    void Set_Beta           (TVector3 const &);
    void Set_P              (TLorentzVector const &);
    void Set_P              (double const &);

    ReactionFragment const & GetReactionFragment(const int & nr) const {
        switch (nr) {
            case 1:
                return p1;
            case 2:
                return p2;
            case 3:
                return p3;
            case 4:
                return p4;
            default:
                throw std::runtime_error("Reaction fragment chosen not existent\n");
        }
    };
protected:
    ReactionFragment p3;
    ReactionFragment p4;
    double p_f_cm;     //Final impulse in COM
    //double p_i;   //Initial impulse in COM
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