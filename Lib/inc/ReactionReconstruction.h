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
    T precision{1.E-8};
public:
    ReactionReconstruction()=delete;
    virtual void SetBeamEnergy(T const &);
    virtual void SetPrecision(T const&);
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
    void SetPrecision(T const & E) override;

    void Set_E              (T const &);
    void Set_Ek             (T const &);
    T Set_E_cm              (T const &);
    T Set_Ek_cm             (T const &);
    void Set_Theta          (T const &, bool);
    void Set_Theta_cm       (T const &);
    void Set_Theta_Phi      (T const &, T const &, bool);
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

public:
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

//IMPLEMENTATION///////////////////////////////////////////////////////////////////////////

//Base cass////////////////////////////////////////////////////////////////////////////////
template<typename T>
ReactionReconstruction<T>::ReactionReconstruction( ReactionFragment::FragmentSettings const & p1,
                                                   ReactionFragment::FragmentSettings const & p2):
        p1(p1),
        p2(p2),
        P_TOT(TLorentzVector()),
        changed_initial_conditions(true){
}

template<typename T>
void ReactionReconstruction<T>::SetBeamEnergy(const T & Energy) {
    p1.Set_Ek(Energy);
    changed_initial_conditions=true;
    UpdateVectors();
}

template<typename T>
void ReactionReconstruction<T>::SetPrecision(const T & precision){
    this->precision = precision;
}

template<typename T>
void ReactionReconstruction<T>::UpdateVectors() {
    P_TOT = p1.Get_P() + p2.Get_P();
    p1.Set_betacm(P_TOT.BoostVector());
    p1.BoostToCm();
    p2.Set_betacm(P_TOT.BoostVector());
    p2.BoostToCm();
}

//2 body cass////////////////////////////////////////////////////////////////////////////////
template<typename T>
ReactionReconstruction2body<T>::ReactionReconstruction2body(ReactionInput2body const & data):
        ReactionReconstruction<T>(data[0], data[1]),
        p3(data[2]),
        p4(data[3]){
    ReactionReconstruction<T>::UpdateVectors();
    ChooseFixed(3);     //Will be computing data based on E3 or Theta3
    ChooseExFixed(3);   //Will be computing Ex4
    if (!ReactionReconstruction2body::CheckConsistency())
        throw std::runtime_error("Data in 2 body reaction initialization not consistent\n");
}

template<typename T>
std::string ReactionReconstruction2body<T>::Get_Name() const{
    return  this->p1.Get_Name()+
            "("+
            this->p2.Get_Name()+
            ","+
            p3.Get_Name()+
            ")"+
            p4.Get_Name()+
            "@"+
            std::to_string(this->p1.Get_Ek()/UNITS::MeV);
}

template<typename T>
void ReactionReconstruction2body<T>::ChooseFixed(const int & nr) {
    if ( nr == 3 ) {
        p3.Set_Fixed(true);
        p4.Set_Fixed(false);
        return;
    }
    if ( nr == 4 ) {
        p3.Set_Fixed(false);
        p4.Set_Fixed(true);
    } else
        throw std::runtime_error("Setting fixed reaction fragment not consistent\n");
}

template<typename T>
void ReactionReconstruction2body<T>::ChooseExFixed(const int & nr) {
    if ( nr == 3 ) {
        p3.Set_ExFixed(true);
        p4.Set_ExFixed(false);
        return;
    }
    if ( nr == 4 ) {
        p3.Set_ExFixed(false);
        p4.Set_ExFixed(true);
    } else
        throw std::runtime_error("Setting Ex fixed reaction fragment not consistent\n");
}

template<typename T>
ReactionFragment& ReactionReconstruction2body<T>::GetFreeFragment() {
    if ( !p3.Is_Fixed() )
        return p3;
    else
    if ( !p4.Is_Fixed() )
        return p4;
    else
        throw std::runtime_error("Inconsistent Fixed reaction fragment\n");
}

template<typename T>
ReactionFragment& ReactionReconstruction2body<T>::GetFixedFragment() {
    if ( p3.Is_Fixed() )
        return p3;
    else
    if ( p4.Is_Fixed() )
        return p4;
    else
        throw std::runtime_error("Inconsistent Fixed reaction fragment\n");
}

template<typename T>
ReactionFragment& ReactionReconstruction2body<T>::GetExFreeFragment() {
    if ( !p3.Is_ExFixed() )
        return p3;
    else
    if ( !p4.Is_ExFixed() )
        return p4;
    else
        throw std::runtime_error("Inconsistent Fixed ex reaction fragment\n");
}

template<typename T>
ReactionFragment& ReactionReconstruction2body<T>::GetExFixedFragment() {
    if ( p3.Is_ExFixed() )
        return p3;
    else
    if ( p4.Is_ExFixed() )
        return p4;
    else
        throw std::runtime_error("Inconsistent Fixed ex reaction fragment\n");
}

template<typename T>
T ReactionReconstruction2body<T>::Get_ThetaMax() {//Checked
    auto & fixed = GetFixedFragment();
    auto & free = GetFreeFragment();
    T val1 = 2.*fixed.Get_M()*(this->p1.Get_E()+this->p2.Get_M());
    T val2 = free.Get_M2()- this->p1.Get_M2() - this->p2.Get_M2() - fixed.Get_M2() -2.*this->p1.Get_E()* this->p2.Get_M();
    T sqroot = val1*val1-val2*val2;
    if (sqroot<0){
        return UNITS::CONSTANTS::pi;
    }else{
        sqroot = sqrtl(sqroot)/(2.*fixed.Get_M()*this->p1.Get_P().Vect().Mag());
        return acosl(sqroot);
    }
}

template<typename T>
void ReactionReconstruction2body<T>::Set_E(const T & E) {//Checked!
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();

    if (E< fixed.Get_M())
        throw std::runtime_error("Relativistic energy is less than mass!\n");

    T eFree = this->p1.Get_E() + this->p2.Get_M() - E;

    fixed.Set_E_Theta(  E,
                        acosl(
                                (this->p2.Get_M2() + free.Get_M2() - this->p1.Get_M2() - fixed.Get_M2() -2* this->p2.Get_M()*(eFree) + 2 * this->p1.Get_E() * E) /
                                sqrt(4*(this->p1.Get_E()*this->p1.Get_E()-this->p1.Get_M2())*(E*E-fixed.Get_M2()))
                        ));

    free.Set_E_Theta_Phi(eFree,
                         acosl(
                                 (this->p2.Get_M2()+fixed.Get_M2()-this->p1.Get_M2()-free.Get_M2()-2* this->p2.Get_M()*(E)+ 2 * this->p1.Get_E() * eFree) /
                                 sqrt(4*(this->p1.Get_E()*this->p1.Get_E()-this->p1.Get_M2())*(eFree * eFree - free.Get_M2()))),
                           UNITS::CONSTANTS::pi+fixed.Get_P().Vect().Phi()
    );
    CheckVectors();
    UpdateVectors(true);
}

template<typename T>
void ReactionReconstruction2body<T>::Set_Ek(const T & Ek) {//Checked
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();
    Set_E(Ek+fixed.Get_M());
}

//Returns Ex value
template<typename T>
T ReactionReconstruction2body<T>::Set_E_cm(const T & E) {
    T ex;
    T s = powl(this->p1.Get_E_cm() + this->p2.Get_E_cm(), 2);
    T sqrtS = this->p1.Get_E_cm() + this->p2.Get_E_cm();
    auto& fixed = GetFixedFragment();
    auto& free = GetFreeFragment();
    if(GetExFixedFragment() == fixed){//Ex4 from E3 -> works
        ex = sqrtl(s+fixed.Get_M2()- 2. * sqrtS * E) - GetExFreeFragment().Get_M();
    }else{//Ex3 from E3 ->does not work with Ek
        ex = sqrtl(-s+free.Get_M2()+ 2. * sqrtS * E) - GetExFreeFragment().Get_M();
    }
    GetExFreeFragment().Set_Ex(ex);
    fixed.Set_E_cm(E);
    free.Set_E_cm((s+free.Get_M2()-fixed.Get_M2())/(2. * sqrtS));
    UpdateVectors(true);
    return ex;
}

template<typename T>
T ReactionReconstruction2body<T>::Set_Ek_cm(const T &Ek) {
    if(GetFixedFragment()==GetExFixedFragment()) {
        return Set_E_cm(Ek + GetFixedFragment().Get_M());
    }else{
        auto& fixed = GetFixedFragment();
        auto& free = GetFreeFragment();
        T sqrtS = this->p1.Get_E_cm() + this->p2.Get_E_cm();
        T s = powl(this->p1.Get_E_cm() + this->p2.Get_E_cm(), 2);
        T ex = -sqrtl(free.Get_M2()+ 2. * Ek * sqrtS) + sqrtS - fixed.Get_M();
        GetExFreeFragment().Set_Ex(ex);
        fixed.Set_E_cm(Ek+fixed.Get_M());
        free.Set_E_cm((s+free.Get_M2()-fixed.Get_M2())/(2. * sqrtS));
        UpdateVectors(true);
        return ex;
    }
}


template<typename T>
void ReactionReconstruction2body<T>::Set_Theta(const T & Theta, bool choseMaxSolution) {//Checked
    double thetaMax = Get_ThetaMax();
    if (thetaMax==UNITS::CONSTANTS::pi) {
        if (Theta < UNITS::CONSTANTS::pi / 2) {
            choseMaxSolution = true;
        }else {
            choseMaxSolution = false;
        }
    }else{
        if (Theta > thetaMax)
            throw std::runtime_error("Theta is exceding theta max\n");
    }

    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();

    T m12 = this->p1.Get_M2();
    T m22 = this->p2.Get_M2();
    T m2 = this->p2.Get_M();
    T m32 = fixed.Get_M2();
    //T m3 = fixed.getM();
    T m42 = free.Get_M2();
    T e1 = this->p1.Get_E();

    T cosine = cos(Theta);
    T p1Mag2 = this->p1.Get_P().Vect().Mag2();
    T sqroot = cosine * cosine * p1Mag2 *
               (     m12*m12+m22*m22+m32*m32+m42*m42+
                      -2.*(m22*m32+m22*m42+m32*m42+m12*m42-m12*m22-m12*m32)
                      + 4. * e1 * e1 * (m22 - m32)
                      + 4. * e1 * m2 * (m12 + m22 - m32 - m42)
                      + 4. * m32 * p1Mag2 * cosine * cosine
                );

    if (sqroot<0)
        throw std::runtime_error("No solution exists: no angle possible\n");
    sqroot = sqrtl(sqroot);

    T numerator = (e1 + m2) * (m12 + m22 + m32 - m42 + 2 * m2 * e1);
    T denominator = 2. * (powl(e1 + m2, 2) - p1Mag2 * cosine * cosine);
    T e3 = 0;

    if (choseMaxSolution)
        e3 = (numerator+sqroot)/denominator;
    else
        e3 = (numerator-sqroot)/denominator;
    Set_E(e3);
}

template<typename T>
void ReactionReconstruction2body<T>::Set_Theta_cm(const T & Theta) {//Not yet correct
    auto &fixed = GetFixedFragment();
    auto &free = GetFreeFragment();
    T s = powl(this->p1.Get_E_cm() + this->p2.Get_E_cm(), 2);
    T pf2 = 0.25/s *(s-powl(p3.Get_M()-p4.Get_M(), 2))*(s-powl(p3.Get_M()+p4.Get_M(), 2));
    T enFix = sqrtl(pf2 + fixed.Get_M2());
    T enFree = sqrtl(pf2 + free.Get_M2());
    fixed.Set_E_Theta_cm(enFix, Theta);
    free.Set_E_Theta_cm(enFree, UNITS::CONSTANTS::pi + Theta);
    UpdateVectors(false);
}

template<typename T>
void ReactionReconstruction2body<T>::Set_Theta_Phi(const T & Theta,
                                                   const T & Phi,
                                                   bool choseMaxSolution) {
    GetFixedFragment().Set_Theta_Phi(Theta, Phi);
    Set_Theta(Theta, choseMaxSolution);
}

template<typename T>
void ReactionReconstruction2body<T>::Set_Theta_Phi_cm(const T & Theta, const T & Phi) {
    GetFixedFragment().Set_Theta_Phi_cm(Theta, Phi);
    GetFixedFragment().Set_Theta_cm(Theta);
}

//Returns Ex value
template<typename T>
T ReactionReconstruction2body<T>::Set_E_Theta(const T & E,
                                              const T & Theta,
                                              const bool& angle_from_fixed) {
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();
    bool chose_max = true;

    GetExFreeFragment().Set_Ex(0);//Should not be needed
    if(angle_from_fixed) {
        if (fixed.Is_ExFixed()) {//Compute Ex4 from E3 and Theta3 -> correct
            T costheta = cosl(Theta);
            T p1mag = this->p1.Get_P().Vect().Mag();
            T sqrt1 = E*E - fixed.Get_M2();
            if (sqrt1 < 0) throw std::runtime_error("negative first square root in ex computation!\n");
            T arg1 = this->p1.Get_M2() + this->p2.Get_M2() + fixed.Get_M2()-2.*E*this->p2.Get_M()+2.*this->p1.Get_E()*(this->p2.Get_M()-E);
            T sqrt2 = arg1 + 2.*costheta*p1mag* sqrtl(sqrt1);
            if (sqrt2 < 0) throw std::runtime_error("negative second square root in ex computation!\n");
            T ex = sqrtl(sqrt2)-GetExFreeFragment().Get_M();
            GetExFreeFragment().Set_Ex(ex);
            Set_E(E);
            if (abs(Theta - fixed.Get_P().Theta())> this->precision)
                throw std::runtime_error ("Mismatch between computed and provided thetas\n");
            UpdateVectors(true);
            return ex;
        } else {//Compute Ex3 from E3 and Theta3 -> not correct with Ek
            T costheta2 = powl(cosl(Theta), 2);
            T p1_mag2 = this->p1.Get_P().Vect().Mag2();
            T cos2_p12 = costheta2*p1_mag2;
            T sqrt1 = cos2_p12*(this->p1.Get_M2()+(E-this->p2.Get_M())*(E-2.*this->p1.Get_E()-this->p2.Get_M())-free.Get_M2()+cos2_p12);
            if (sqrt1 < 0) throw std::runtime_error("negative first square root in ex computation!\n");
            T arg1 = 2.*this->p1.Get_E()*(E-this->p2.Get_M())+2.*E*this->p2.Get_M()-this->p1.Get_M2()-this->p2.Get_M2()+free.Get_M2();
            if (!chose_max){
                T sqrt2 = arg1 - 2.*(cos2_p12+sqrtl(sqrt1));
                std::cout << "sqrt 2 " << sqrt2 << std::endl;
                if (sqrt2<0) {
                    std::cout << "switching to max solution\n";
                    chose_max = true;
                }else {
                    T ex = sqrtl(sqrt2)-GetExFreeFragment().Get_M();
                    GetExFreeFragment().Set_Ex(ex);
                    Set_E(E);
                    if (abs(Theta - fixed.Get_P().Theta())> this->precision)
                        throw std::runtime_error ("Mismatch between computed and provided thetas\n");
                    UpdateVectors(true);
                    return ex;
                }
            }
            if (chose_max){
                T sqrt2 = arg1 - 2.*(cos2_p12-sqrtl(sqrt1));
                if (sqrt2 < 0) throw std::runtime_error("negative second square root in ex computation!\n");
                T ex = sqrtl(sqrt2)-GetExFreeFragment().Get_M();
                GetExFreeFragment().Set_Ex(ex);
                Set_E(E);
                if (abs(Theta - fixed.Get_P().Theta())> this->precision)
                    throw std::runtime_error ("Mismatch between computed and provided thetas\n");
                UpdateVectors(true);
                return ex;
            }
        }
    }else{
        if (fixed.Is_ExFixed()) {//Compute Ex4 from E3 and Theta4 -> correct
            T costheta2 = powl(cosl(Theta), 2);
            T p1_mag2 = this->p1.Get_P().Vect().Mag2();
            T cos2_p12 = costheta2*p1_mag2;
            T sqrt1 = cos2_p12*(this->p1.Get_M2()-fixed.Get_M2()-this->p1.Get_E()*this->p1.Get_E()+E*E+cos2_p12);
            if (sqrt1 < 0) throw std::runtime_error("negative first square root in ex computation!\n");
            T arg1 = -this->p1.Get_M2()+this->p2.Get_M2()+this->p3.Get_M2()+2.*(this->p1.Get_E()-E)*(this->p1.Get_E()+this->p2.Get_M());
            if (!chose_max){
                T sqrt2 = arg1 - 2.*(cos2_p12+sqrtl(sqrt1));
                if (sqrt2<0) {
                    std::cout << "switching to max solution\n";
                    chose_max = true;
                }else {
                    T ex = sqrtl(sqrt2)-GetExFreeFragment().Get_M();
                    GetExFreeFragment().Set_Ex(ex);
                    Set_E(E);
                    if (abs(Theta - free.Get_P().Theta())> this->precision)
                        throw std::runtime_error ("Mismatch between computed and provided thetas\n");
                    UpdateVectors(true);
                    return ex;
                }
            }
            if (chose_max){
                T sqrt2 = arg1 - 2.*(cos2_p12-sqrtl(sqrt1));
                if (sqrt2 < 0) throw std::runtime_error("negative second square root in ex computation!\n");
                T ex = sqrtl(sqrt2)-GetExFreeFragment().Get_M();
                GetExFreeFragment().Set_Ex(ex);
                Set_E(E);
                if (abs(Theta - free.Get_P().Theta())> this->precision)
                    throw std::runtime_error ("Mismatch between computed and provided thetas\n");
                UpdateVectors(true);
                return ex;
            }
        }else{//Compute Ex3 from E3 and Theta4 -> Not correct with Ek
            T costheta = cosl(Theta);
            T p1mag = this->p1.Get_P().Vect().Mag();
            T sqrt1 = powl(this->p1.Get_E()+this->p2.Get_M()-E, 2) - free.Get_M2();
            if (sqrt1 < 0) throw std::runtime_error("negative first square root in ex computation!\n");
            T arg1 = this->p1.Get_M2() - this->p2.Get_M2() + free.Get_M2()-2.*powl(this->p1.Get_E(), 2)+2.*this->p1.Get_E()*(E-this->p2.Get_M());
            T sqrt2 = arg1 + 2.*costheta*p1mag* sqrtl(sqrt1);
            if (sqrt2 < 0) throw std::runtime_error("negative second square root in ex computation!\n");
            T ex = sqrtl(sqrt2)-GetExFreeFragment().Get_M();
            GetExFreeFragment().Set_Ex(ex);
            Set_E(E);
            if (abs(Theta - free.Get_P().Theta())> this->precision)
                throw std::runtime_error ("Mismatch between computed and provided thetas\n");
            UpdateVectors(true);
            return ex;

        }
    }
    return -1;
}

template<typename T>
T ReactionReconstruction2body<T>::Set_Ek_Theta(const T & Ek,
                                               const T & Theta,
                                               const bool& angle_from_fixed) {
    //Returns Ex value
    if (GetExFixedFragment()==GetFixedFragment())
        return Set_E_Theta(Ek+GetFixedFragment().Get_M(), Theta, angle_from_fixed);
    else
        throw
                std::runtime_error("feature not implemented\n");
    //è qui l'errore!!!!
}

template<typename T>
T ReactionReconstruction2body<T>::Set_E_Theta_cm(const T & E, const T & Theta) {
    Set_Theta_cm(Theta);
    return Set_E_cm(E);
}

template<typename T>
T ReactionReconstruction2body<T>::Set_Ek_Theta_cm(const T & E, const T & Theta) {
    Set_Theta_cm(Theta);
    return Set_Ek_cm(E);
}

template<typename T>
void ReactionReconstruction2body<T>::Set_Beta(const TVector3 & Beta) {
    GetFixedFragment().Set_Theta_Phi(Beta.Theta(), Beta.Phi());
    Set_E(GetFixedFragment().Get_M()/sqrtl(1-Beta.Mag2()));
}

template<typename T>
void ReactionReconstruction2body<T>::Set_P(const T & P) {
    Set_E(sqrtl(GetFixedFragment().Get_M2()+P*P));
}

template<typename T>
void ReactionReconstruction2body<T>::CheckVectors(){
    TLorentzVector sum = this->p1.Get_P()+this->p2.Get_P();
    TLorentzVector discr = sum-p3.Get_P()-p4.Get_P();

    if (    abs(discr.E())  > this->precision ||
            (abs(discr.Px())> this->precision )||
            (abs(discr.Px())> this->precision )||
            (abs(discr.Px())> this->precision ))
        throw std::runtime_error("Not enough precision in vectors\n");
}

template<typename T>
bool ReactionReconstruction2body<T>::CheckConsistency() const {
    if(this->p1.Get_A()+this->p2.Get_A()!=p3.Get_A()+p4.Get_A()) {
        std::cerr << "Total A does not match\n";
        return false;
    }
    if(this->p1.Get_Z() + this->p2.Get_Z() != p3.Get_Z() + p4.Get_Z()) {
        std::cerr << "Total Z does not match\n";
        return false;
    }
    return p3.Is_Fixed() != p4.Is_Fixed();
}

template<typename T>
void ReactionReconstruction2body<T>::SetBeamEnergy(const T &E) {
    ReactionReconstruction<T>::SetBeamEnergy(E);
}

template<typename T>
void ReactionReconstruction2body<T>::SetPrecision(const T & precision) {
    ReactionReconstruction<T>::SetPrecision(precision);
}

template<typename T>
T ReactionReconstruction2body<T>::Set_E_Theta_Phi(const T & E,
                                                  const T & Theta,
                                                  const T & Phi) {
    GetFixedFragment().Set_Theta_Phi(Theta, Phi);
    return Set_E_Theta(E, Theta);
}

template<typename T>
T ReactionReconstruction2body<T>::Set_Ek_Theta_Phi(const T & E,
                                                   const T & Theta,
                                                   const T & Phi) {
    GetFixedFragment().Set_Theta_Phi(Theta, Phi);
    return Set_Ek_Theta(E, Theta);
}

template<typename T>
void ReactionReconstruction2body<T>::UpdateVectors(const bool& lab) {
    if (this->changed_initial_conditions) {
        ReactionReconstruction<T>::UpdateVectors();
        p3.Set_betacm(this->p1.Get_betacm());
        p4.Set_betacm(this->p1.Get_betacm());
    }
    if (lab){
        p3.BoostToCm();
        p4.BoostToCm();
    }else{
        p3.BoostToLab();
        p4.BoostToLab();
    }
    this->changed_initial_conditions = false;
}

//3 body cass////////////////////////////////////////////////////////////////////////////////
template<typename T>
ReactionReconstruction3body<T>::ReactionReconstruction3body(std::array<ReactionFragment::FragmentSettings, 4> const & data,
                                                            ReactionFragment::FragmentSettings const & p5,
                                                            ReactionFragment::FragmentSettings const & p6):
        ReactionReconstruction2body<T>(data),
        p5(p5),
        p6(p6){
    if (!ReactionReconstruction3body::CheckConsistency())
        throw std::runtime_error("Data in 3 body reaction initialization not consistent\n");
}

template<typename T>
std::string ReactionReconstruction3body<T>::Get_Name() const{
    return  ReactionReconstruction2body<T>::Get_Name()+
            "*->"+
            p5.Get_Name();
}

template<typename T>
bool ReactionReconstruction3body<T>::CheckConsistency() const {
    if(!ReactionReconstruction2body<T>::CheckConsistency())
        return false;
    if(this->p4.Get_A() != p5.Get_A()+p6.Get_A())
        return false;
    else
        return this->p4.Get_Z() != p5.Get_Z()+p6.Get_Z();
}

