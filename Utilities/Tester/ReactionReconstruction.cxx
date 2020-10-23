#include <ReactionReconstruction.h>

//#define sqrt sqrtl

//Base cass////////////////////////////////////////////////////////////////////////////////
ReactionReconstruction::ReactionReconstruction( ReactionFragment::FragmentSettings const & p1,
                                                ReactionFragment::FragmentSettings const & p2):
        p1(p1),
        p2(p2),
        P_TOT(TLorentzVector()),
        changed_initial_conditions(true){
    //UpdateVectors();
}

void ReactionReconstruction::SetBeamEnergy(const long double & Energy) {
    p1.Set_Ek(Energy);
    changed_initial_conditions = true;
    UpdateVectors();
}

void ReactionReconstruction::UpdateVectors() {
    if (!changed_initial_conditions)
        return;
    P_TOT = p1.Get_P() + p2.Get_P();

//    long double debug1 = P_TOT*P_TOT;
//    std::cout <<"prod: " <<debug1 << std::endl;
//    long double debug = p1.Get_M2()+p2.Get_M2()+2*p1.Get_E()*p2.Get_M();
//    std::cout << "comp : " <<debug << std::endl;
    p1.Set_Invariant(P_TOT*P_TOT);
///    long double debug2 = p1.Get_Invariant();
    p2.Set_Invariant(p1.Get_Invariant());
    p1.Set_betacm(P_TOT.BoostVector()); //Not necessary
    p2.Set_betacm(P_TOT.BoostVector()); //Not necessary
    changed_initial_conditions = false;
}

//2 body cass////////////////////////////////////////////////////////////////////////////////
ReactionReconstruction2body::ReactionReconstruction2body(ReactionInput2body const & data):
        ReactionReconstruction(data[0], data[1]),
        p3(data[2]),
        p4(data[3]){
    ReactionReconstruction::UpdateVectors();
    ChooseFixed(3);     //Will be computing data based on E3 or Theta3
    ChooseExFixed(3);   //Will be computing Ex4
    if (!ReactionReconstruction2body::CheckConsistency())
        throw std::runtime_error("Data in 2 body reaction initialization not consistent\n");
}

std::string ReactionReconstruction2body::Get_Name() const{
    return  p1.Get_Name()+
            "@"+
            std::to_string(p1.Get_Ek()/UNITS::MeV)+
            "("+
            p2.Get_Name()+
            ","+
            p3.Get_Name()+
            ")";
}

void ReactionReconstruction2body::ChooseFixed(const int & nr) {
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

void ReactionReconstruction2body::ChooseExFixed(const int & nr) {
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

ReactionFragment& ReactionReconstruction2body::GetFreeFragment() {
    if ( !p3.Is_Fixed() )
        return p3;
    else
        if ( !p4.Is_Fixed() )
            return p4;
        else
            throw std::runtime_error("Inconsistent Fixed reaction fragment\n");
}

ReactionFragment& ReactionReconstruction2body::GetFixedFragment() {
    if ( p3.Is_Fixed() )
        return p3;
    else
        if ( p4.Is_Fixed() )
            return p4;
        else
            throw std::runtime_error("Inconsistent Fixed reaction fragment\n");
}

ReactionFragment& ReactionReconstruction2body::GetExFreeFragment() {
    if ( !p3.Is_ExFixed() )
        return p3;
    else
    if ( !p4.Is_ExFixed() )
        return p4;
    else
        throw std::runtime_error("Inconsistent Fixed ex reaction fragment\n");
}

ReactionFragment& ReactionReconstruction2body::GetExFixedFragment() {
    if ( p3.Is_ExFixed() )
        return p3;
    else
    if ( p4.Is_ExFixed() )
        return p4;
    else
        throw std::runtime_error("Inconsistent Fixed ex reaction fragment\n");
}

long double ReactionReconstruction2body::Get_ThetaMax() {//Checked
    auto & fixed = GetFixedFragment();
    auto & free = GetFreeFragment();
    long double val1 = 2.*fixed.Get_M()*(p1.Get_E()+p2.Get_M());
    long double val2 = free.Get_M2()- p1.Get_M2() - p2.Get_M2() - fixed.Get_M2() -2.*p1.Get_E()* p2.Get_M();
    long double sqroot = val1*val1-val2*val2;
    if (sqroot<0){
        return UNITS::CONSTANTS::pi;
    }else{
        sqroot = sqrtl(sqroot)/(2.*fixed.Get_M()*p1.Get_P().Vect().Mag());
        return acosl(sqroot);
    }
}

void ReactionReconstruction2body::Set_E(const long double & E) {//Checked!
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();

    if (E< fixed.Get_M())
        throw std::runtime_error("Relativistic energy is less than mass!\n");

    std::cout << "Mass before: " << fixed.Get_M() << std::endl;
    std::cout << "Energy before: " << E << std::endl;
    long double E_free = p1.Get_E()+p2.Get_M()-E;

    fixed.Set_E_Theta(  E,
                         acosl(
                            (p2.Get_M2()+free.Get_M2()-p1.Get_M2()-fixed.Get_M2()-2* p2.Get_M()*(E_free)+2*p1.Get_E()*E)/
                                            sqrt(4*(p1.Get_E()*p1.Get_E()-p1.Get_M2())*(E*E-fixed.Get_M2()))
                         ));

    free.Set_E_Theta_Phi(  E_free,
                           acosl(
                           (p2.Get_M2()+fixed.Get_M2()-p1.Get_M2()-free.Get_M2()-2* p2.Get_M()*(E)+2*p1.Get_E()*E_free)/
                           sqrt(4*(p1.Get_E()*p1.Get_E()-p1.Get_M2())*(E_free*E_free-free.Get_M2()))),
                           UNITS::CONSTANTS::pi+fixed.Get_P().Vect().Phi()
                           );

    std::cout << "fixed phi "<< fixed.Get_P().Vect().X() <<std::endl;
    std::cout << "free phi "<< free.Get_P().Vect().X() <<std::endl;
    //UpdateVectors();
    CheckVectors();
}

void ReactionReconstruction2body::Set_Ek(const long double & Ek) {//Checked
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();

    long double debug = 0;
    debug = p1.Get_M()/UNITS::MeV ;
    debug = p2.Get_M()/UNITS::MeV ;
    debug = p3.Get_M()/UNITS::MeV ;
    debug = p4.Get_M()/UNITS::MeV ;
    Set_E(Ek+fixed.Get_M());
}

//Returns Ex value
long double ReactionReconstruction2body::Set_E_cm(const long double & E) {
    long double Ex = 0;

    if(p3.Is_Fixed() && p3.Is_ExFixed()) {//Compute Ex4 from E3
        Ex = sqrtl(-2 * sqrtl(p1.Get_Invariant()) * E + p1.Get_Invariant() + p3.Get_M2()) - p4.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    if(!p3.Is_Fixed() && p3.Is_ExFixed()) {//Compute Ex4 from E4
        Ex = sqrtl(2 * sqrtl(p1.Get_Invariant()) * E - p1.Get_Invariant() + p3.Get_M2()) - p4.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    if(p3.Is_Fixed() && !p3.Is_ExFixed()) {//Compute Ex3 from E3
        Ex = sqrtl(2 * sqrtl(p1.Get_Invariant()) * E - p1.Get_Invariant() + p4.Get_M2()) - p3.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    if(!p3.Is_Fixed() && !p3.Is_ExFixed()) {//Compute Ex3 from E4
        Ex = sqrtl(-2 * sqrtl(p1.Get_Invariant()) * E + p1.Get_Invariant() + p4.Get_M2()) - p3.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    throw
        std::runtime_error("Something wrong in fixed and ex-fixed data\n");
}

void ReactionReconstruction2body::Set_Theta(const long double & Theta, const bool & chose_max_solution=true) {//Checked
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();
    long double m12 = p1.Get_M2();
    long double m22 = p2.Get_M2();
    long double m2 = p2.Get_M();
    long double m32 = fixed.Get_M2();
    //long double m3 = fixed.Get_M();
    long double m42 = free.Get_M2();
    long double E1 = p1.Get_E();

    long double cosine = cos(Theta);
    long double p1_mag2 = p1.Get_P().Vect().Mag2();
    //free.Set_Invariant(1);
    long double sqroot =  cosine*cosine*p1_mag2 *
                            (     m12*m12+m22*m22+m32*m32+m42*m42+
                             -2.*(m22*m32+m22*m42+m32*m42+m12*m42-m12*m22-m12*m32)
                             +4.*E1*E1*(m22-m32)
                             +4.*E1*m2*(m12+m22-m32-m42)
                             +4.*m32*p1_mag2*cosine*cosine
                            );

    if (sqroot<0)
        throw std::runtime_error("No solution exists: no angle possible\n");

    sqroot = sqrtl(sqroot);

    long double numerator = (E1+m2)*(m12+m22+m32-m42+2*m2*E1);
    long double denominator = 2 * ( pow(E1+m2, 2)-p1_mag2*cosine*cosine);
    long double e3 = 0;

    if (chose_max_solution)
        e3 = (numerator+sqroot)/denominator;
    else
        e3 = (numerator-sqroot)/denominator;

    Set_E(e3);
}

void ReactionReconstruction2body::Set_Theta_cm(const long double & Theta) {
    GetFixedFragment().Set_Theta_cm(Theta);
    UpdateVectors();
}

void ReactionReconstruction2body::Set_Theta_Phi(const long double & Theta, const long double & Phi, bool const & choose_max) {
    GetFixedFragment().Set_Theta_Phi(Theta, Phi);
    Set_Theta(Theta, choose_max);
}

void ReactionReconstruction2body::Set_Theta_Phi_cm(const long double & Theta, const long double & Phi) {
    GetFixedFragment().Set_Theta_Phi_cm(Theta, Phi);
    UpdateVectors();
}

//Returns Ex value
long double ReactionReconstruction2body::Set_E_Theta(const long double & E, const long double & Theta, bool chose_max = false, const bool& angle_from_fixed= true) {
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();
    long double Ex = 0;

    ///SET FREE EX TO ZERO!!!!!!!
    GetExFreeFragment().Set_Ex(0);
    if(angle_from_fixed) {
        if (fixed.Is_ExFixed()) {//Compute Ex4 from E3 and Theta3
            long double costheta = cosl(Theta);
            long double p1mag = p1.Get_P().Vect().Mag();
            long double sqrt1 = E*E - fixed.Get_M2();
            if (sqrt1 < 0) throw std::runtime_error("negative first square root in ex computation!\n");
            long double arg1 = p1.Get_M2() + p2.Get_M2() + fixed.Get_M2()-2.*E*p2.Get_M()+2.*p1.Get_E()*(p2.Get_M()-E);
            long double sqrt2 = arg1 + 2.*costheta*p1mag* sqrtl(sqrt1);
            if (sqrt2 < 0) throw std::runtime_error("negative second square root in ex computation!\n");
            long double p = sqrtl(E*E-fixed.Get_M2());
            TLorentzVector a(p*sinl(Theta), 0., p*cosl(Theta), E);
            std::cout << "should be : " << (p1.Get_P()+p2.Get_P()-a).Mag() - GetExFreeFragment().Get_M() << std::endl;
            return sqrtl(sqrt2)-GetExFreeFragment().Get_M();

        } else {//Compute Ex3 from E3 and Theta3
            long double costheta2 = pow(cosl(Theta), 2);
            long double p1_mag2 = p1.Get_P().Vect().Mag2();
            long double cos2_p12 = costheta2*p1_mag2;
            long double sqrt1 = cos2_p12*(p1.Get_M2()+(E-p2.Get_M())*(E-2.*p1.Get_E()-p2.Get_M())-free.Get_M2()+cos2_p12);
            std::cout << "sqrt 1 " << sqrt1 << std::endl;
            if (sqrt1 < 0) throw std::runtime_error("negative first square root in ex computation!\n");
            long double arg1 = 2.*p1.Get_E()*(E-p2.Get_M())+2.*E*p2.Get_M()-p1.Get_M2()-p2.Get_M2()+free.Get_M2();
            std::cout << "arg 1 " << arg1 << std::endl;
            if (!chose_max){
                long double sqrt2 = arg1 - 2.*(cos2_p12+sqrtl(sqrt1));
                std::cout << "sqrt 2 " << sqrt2 << std::endl;
                if (sqrt2<0) {
                    std::cout << "switching to max solution\n";
                    chose_max = true;
                }else {
                    std::cout << "res sqrt 2 " << sqrtl(sqrt2) << std::endl;
                    std::cout << "res ex " << sqrtl(sqrt2) -GetExFreeFragment().Get_M()<< std::endl;
                    return sqrtl(sqrt2)-GetExFreeFragment().Get_M();
                }
            }
            if (chose_max){
                long double sqrt2 = arg1 - 2.*(cos2_p12-sqrtl(sqrt1));
                if (sqrt2 < 0) throw std::runtime_error("negative second square root in ex computation!\n");
                return sqrtl(sqrt2)-GetExFreeFragment().Get_M();
            }
        }
    }else{
        if (fixed.Is_ExFixed()) {//Compute Ex4 from E3 and Theta4

        }else{//Compute Ex3 from E3 and Theta4

        }
    }
}

long double ReactionReconstruction2body::Set_Ek_Theta(const long double & Ek, const long double & Theta, bool chose_max = false, const bool& angle_from_fixed= true) {
    return Set_E_Theta(Ek+GetFixedFragment().Get_M(), Theta, chose_max, angle_from_fixed);
}
//Returns Ex value
long double ReactionReconstruction2body::Set_E_Theta_cm(const long double & E, const long double & Theta) {
    Set_Theta_cm(Theta);
    return Set_E_cm(E);
}

void ReactionReconstruction2body::Set_Dir(const TVector3 & Dir) {
    GetFreeFragment().Set_Dir(Dir);
    UpdateVectors();
}

void ReactionReconstruction2body::Set_Beta(const TVector3 & Beta) {
    GetFreeFragment().Set_Beta(Beta);
    UpdateVectors();
}

void ReactionReconstruction2body::Set_P(const TLorentzVector & P) {
    GetFreeFragment().Set_P(P);
    UpdateVectors();
}

void ReactionReconstruction2body::Set_P(const long double & P) {
    Set_E(sqrtl(GetFixedFragment().Get_M2()+P*P));
}

void ReactionReconstruction2body::UpdateVectors() {
    if (changed_initial_conditions) {
        ReactionReconstruction::UpdateVectors();
        p3.Set_betacm(P_TOT.BoostVector());
        p4.Set_betacm(P_TOT.BoostVector());
    }
    GetFreeFragment().Set_P(p1.Get_P()+p2.Get_P()-GetFixedFragment().Get_P());
}

void ReactionReconstruction2body::CheckVectors(){
    TLorentzVector sum = p1.Get_P()+p2.Get_P();
    TLorentzVector discr = sum-p3.Get_P()-p4.Get_P();

//    std::cout << "discr :" <<  abs(discr.E()) << std::endl;
//    std::cout << "discr :" <<  abs(discr.Px()) << std::endl;
//    std::cout << "discr :" <<  abs(discr.Py()) << std::endl;
//    std::cout << "discr :" <<  abs(discr.Pz()) << std::endl;
//    std::cout << "sum :" <<  abs(sum.E()) << std::endl;
//    std::cout << "sum :" <<  abs(sum.Px()) << std::endl;
//    std::cout << "sum :" <<  abs(sum.Py()) << std::endl;
//    std::cout << "sum :" <<  abs(sum.Pz()) << std::endl;

    if (    abs(discr.E())  > precision ||
           (abs(discr.Px())> precision )||
           (abs(discr.Px())> precision )||
           (abs(discr.Px())> precision ))
        throw std::runtime_error("Not enough precision in vectors\n");
}

bool ReactionReconstruction2body::CheckConsistency() const {
    if(p1.Get_A()+p2.Get_A()!=p3.Get_A()+p4.Get_A())
        return false;
    if(p1.Get_Z() + p2.Get_Z() != p3.Get_Z() + p4.Get_Z())
        return false;
    return p3.Is_Fixed() != p4.Is_Fixed();
}

void ReactionReconstruction2body::SetBeamEnergy(const long double &E) {
    ReactionReconstruction::SetBeamEnergy(E);
}

//3 body cass////////////////////////////////////////////////////////////////////////////////
ReactionReconstruction3body::ReactionReconstruction3body(ReactionInput2body const & data,
                                                         ReactionFragment::FragmentSettings const & p5,
                                                         ReactionFragment::FragmentSettings const & p6):
        ReactionReconstruction2body(data),
        p5(p5),
        p6(p6){
    if (!ReactionReconstruction3body::CheckConsistency())
        throw std::runtime_error("Data in 3 body reaction initialization not consistent\n");
}

std::string ReactionReconstruction3body::Get_Name() const{
    return  ReactionReconstruction2body::Get_Name()+
            "*->"+
            p5.Get_Name();
}

bool ReactionReconstruction3body::CheckConsistency() const {
    if(!ReactionReconstruction2body::CheckConsistency())
        return false;
    if(p4.Get_A() != p5.Get_A()+p6.Get_A())
        return false;
    else
        return p4.Get_Z() != p5.Get_Z()+p6.Get_Z();
}
