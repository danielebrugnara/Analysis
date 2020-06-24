#include <ReactionReconstruction.h>

//Base cass////////////////////////////////////////////////////////////////////////////////
ReactionReconstruction::ReactionReconstruction( ReactionFragment::FragmentSettings const & p1,
                                                ReactionFragment::FragmentSettings const & p2):
        p1(p1),
        p2(p2),
        P_TOT(TLorentzVector()),
        changed_initial_conditions(true){
}

void ReactionReconstruction::SetBeamEnergy(const double & Energy) {
    p1.Set_Ek(Energy);
    changed_initial_conditions = true;
    UpdateVectors();
}

void ReactionReconstruction::UpdateVectors() {
    if (!changed_initial_conditions)
        return;
    P_TOT = p1.Get_P() + p2.Get_P();
    p1.Set_Invariant( P_TOT * P_TOT );
    p2.Set_Invariant(p1.Get_Invariant());
    //p1.Set_P_tot(P_CM); //Not necessary
    //p2.Set_P_tot(P_CM); //Not necessary
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

double ReactionReconstruction2body::Get_ThetaMax() {
    p3.Set_Invariant(p2.Get_M2()+p4.Get_M2()-2* p2.Get_M() * p4.Get_E());
    double val1 = 2 * p1.Get_E() * p3.Get_M();
    double val2 = p1.Get_M2()+p3.Get_M2()-p3.Get_Invariant();
    double sqroot = val1*val1 - val2 * val2;
    if (sqroot<0){
        return UNITS::CONSTANTS::pi;
    }else{
        sqroot = sqrt(sqroot);
        return sqroot/(2 * p3.Get_M() * p1.Get_P().Vect().Mag());
    }
}

void ReactionReconstruction2body::Set_E(const double & E) {
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();

    free.Set_Invariant(p2.Get_M2()+fixed.Get_M2()-2*(E)*p2.Get_M());
    fixed.Set_Invariant(2*p1.Get_M2()+2*p2.Get_M2()-free.Get_Invariant()-p1.Get_Invariant());
    fixed.Set_E_Theta(  E,
                        acos(   (fixed.Get_Invariant()-p1.Get_M2()-fixed.Get_M2()+2*E*p1.Get_E())/
                            (2*p1.Get_P().Vect().Mag()*sqrt(E*E-fixed.Get_M2()))));
    UpdateVectors();
}

void ReactionReconstruction2body::Set_Ek(const double & Ek) {
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();

    double E = Ek + fixed.Get_M();
    free.Set_Invariant(p2.Get_M2()+fixed.Get_M2()-2*(E)*p2.Get_M());
    fixed.Set_Invariant(2*p1.Get_M2()+2*p2.Get_M2()-free.Get_Invariant()-p1.Get_Invariant());
    fixed.Set_E_Theta(  E,
                        acos(   (fixed.Get_Invariant()-p1.Get_M2()-fixed.Get_M2()+2*E*p1.Get_E())/
                                (2*p1.Get_P().Vect().Mag()*sqrt(E*E-fixed.Get_M2()))));
    bool debug = fixed.Check_Consistent();
    UpdateVectors();
}

//Returns Ex value
double ReactionReconstruction2body::Set_E_cm(const double & E) {
    double Ex = 0;

    if(p3.Is_Fixed() && p3.Is_ExFixed()) {//Compute Ex4 from E3
        Ex = sqrt(-2 * sqrt(p1.Get_Invariant()) * E + p1.Get_Invariant() + p3.Get_M2()) - p4.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    if(!p3.Is_Fixed() && p3.Is_ExFixed()) {//Compute Ex4 from E4
        Ex = sqrt(2 * sqrt(p1.Get_Invariant()) * E - p1.Get_Invariant() + p3.Get_M2()) - p4.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    if(p3.Is_Fixed() && !p3.Is_ExFixed()) {//Compute Ex3 from E3
        Ex = sqrt(2 * sqrt(p1.Get_Invariant()) * E - p1.Get_Invariant() + p4.Get_M2()) - p3.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    if(!p3.Is_Fixed() && !p3.Is_ExFixed()) {//Compute Ex3 from E4
        Ex = sqrt(-2 * sqrt(p1.Get_Invariant()) * E + p1.Get_Invariant() + p4.Get_M2()) - p3.Get_M_GS();
        p4.Set_Ex(Ex);
        UpdateVectors();
        return Ex;
    }
    throw
        std::runtime_error("Something wrong in fixed and ex-fixed data\n");
}

void ReactionReconstruction2body::Set_Theta(const double & Theta, const bool & chose_max) {
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();

    double cosine = cos(Theta);
    double p1_mag2 = p1.Get_P().Vect().Mag2();
    free.Set_Invariant(1);
    double sqroot =   cosine * cosine * p1_mag2 *
                    ( p1.Get_M2() * p1.Get_M2() -
                        4 * p1.Get_E() * p1.Get_E() * fixed.Get_M2()+
                        fixed.Get_M2() * fixed.Get_M2() +
                        4 * cosine * cosine * fixed.Get_M2() * p1_mag2 +
                        2 * p1.Get_M2() * ( fixed.Get_M2() - fixed.Get_Invariant()) -
                        2 * fixed.Get_M2() * fixed.Get_Invariant() +
                        fixed.Get_Invariant() * fixed.Get_Invariant());
    sqroot = sqrt(sqroot);

    if (! isnan(sqroot))
        return;

    double numerator = 2 * p1.Get_E() * (p1.Get_M2() + fixed.Get_M2() - fixed.Get_Invariant());
    double denominator = 2 * ( p1.Get_E() * p1.Get_E() - cosine * cosine * p1_mag2);
    double e3 = 0;

    if (chose_max)
        e3 = (numerator+sqroot)/denominator;
    else
        e3 = (numerator-sqroot)/denominator;

    Set_E(e3);
}

void ReactionReconstruction2body::Set_Theta_cm(const double & Theta) {
    GetFixedFragment().Set_Theta_cm(Theta);
    UpdateVectors();
}

void ReactionReconstruction2body::Set_Theta_Phi(const double & Theta, const double & Phi, bool const & choose_max) {
    GetFixedFragment().Set_Theta_Phi(Theta, Phi);
    Set_Theta(Theta, choose_max);
}

void ReactionReconstruction2body::Set_Theta_Phi_cm(const double & Theta, const double & Phi) {
    GetFixedFragment().Set_Theta_Phi_cm(Theta, Phi);
    UpdateVectors();
}

//Returns Ex value
double ReactionReconstruction2body::Set_E_Theta(const double & E, const double & Theta, const bool & chose_max = false) {
    auto & fixed    = GetFixedFragment();
    auto & free     = GetFreeFragment();
    double Ex = 0;
    double costheta = cos(Theta);
    double p1_mag = p1.Get_P().Vect().Mag();

    if(fixed.Is_ExFixed()) {//Compute Ex4 from E3 or Ex3 from E4
        if(free.Is_ExFixed())
            throw std::runtime_error("Both Ex are fixed\n");
        fixed.Set_Invariant(   2*p1_mag*sqrt(E*E-fixed.Get_M2())*costheta+
                            p1.Get_M2()+
                            fixed.Get_M2()-
                            2*p1.Get_E()*E);
        //p4.Set_Invariant(2*p1.Get_M2()+p2.Get_M2()-p3.Get_Invariant());
        Ex = sqrt(  2* p2.Get_M()*(  p1.Get_E()+
                                        p2.Get_M()-E)-
                       p2.Get_M2()+
                       fixed.Get_Invariant())
             -free.Get_M_GS();
        free.Set_Ex(Ex);
        fixed.Set_E_Theta( E, Theta);
                        //acos(   (p3.Get_Invariant()-p1.Get_M2()-p3.Get_M2()+2*E*p1.Get_E())/
                        //            (2*p1.Get_P().Vect().Mag()*sqrt(E*E-p3.Get_M2()))));
        UpdateVectors();
        return Ex;
    }else{//Compute Ex3 from E3 or Ex4 from E4
        double efree = p1.Get_E()+p2.Get_E()-E;
        double sqroot = costheta * costheta * p1_mag *p1_mag * ( E*E + p1.Get_M2() +
                                                                2*efree*p2.Get_M() -
                                                                p2.Get_M2() - free.Get_M2() +
                                                                costheta*costheta*p1_mag*p1_mag);
        sqroot = sqrt(sqroot);
        double argm = - p1.Get_M2() - 2*efree*p2.Get_M()+p2.Get_M2()+free.Get_M2()-2*costheta*costheta*p1_mag*p1_mag;

        if(chose_max)
            Ex = sqrt(argm + 2 * sqroot)-fixed.Get_M_GS();
        else
            Ex = sqrt(argm + 2 * sqroot)-fixed.Get_M_GS();
        fixed.Set_Ex(Ex);
        Set_E(E);
        return Ex;
    }
}

//Returns Ex value
double ReactionReconstruction2body::Set_E_Theta_cm(const double & E, const double & Theta) {
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

void ReactionReconstruction2body::Set_P(const double & P) {
    Set_E(sqrt(GetFixedFragment().Get_M2()+P*P));
}

void ReactionReconstruction2body::UpdateVectors() {
    if (changed_initial_conditions) {
        ReactionReconstruction::UpdateVectors();
        p3.Set_P_tot(P_TOT);
        p4.Set_P_tot(P_TOT);
        p3.Set_P_tot_cm(-P_TOT);
        p4.Set_P_tot_cm(-P_TOT);
    }

    //p_f_cm =    (p1.Get_Invariant() - pow(p3.Get_M() - p4.Get_M(), 2) ) /
    //            (2 * p1.Get_Invariant());
    //GetFixedFragment().Set_P(p_f_cm);
    GetFreeFragment().Set_P(p1.Get_P()+p2.Get_P()-GetFixedFragment().Get_P());
}

bool ReactionReconstruction2body::CheckConsistency() const {
    if(p1.Get_A()+p2.Get_A()!=p3.Get_A()+p4.Get_A())
        return false;
    if(p1.Get_Z() + p2.Get_Z() != p3.Get_Z() + p4.Get_Z())
        return false;
    return p3.Is_Fixed() != p4.Is_Fixed();
}

void ReactionReconstruction2body::SetBeamEnergy(const double &E) {
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
