#include <ReactionFragment.h>

#include <iostream>

ReactionFragment::ReactionFragment(FragmentSettings const & data):
        A           (data.A),
        Z           (data.Z),
        Q           (data.Q),
        Spin        (0),
        LifeTime    (-1),
        MassExcess  (-1),
        Parity      (""),
        Name        (""),
        Ex          (data.Ex),
        Lab         (),
        Cm          (),
        Fixed       (true),
        ExFixed     (true),
        precision   (1E-8)
{
    GetFromData(A, Z);
    UpdateMass();
    long double energy = data.Ek + M;
    Lab.P = TLorentzVector( 0.,
                            0.,
                            sqrtl(energy*energy-M2),
                            energy);
}

void ReactionFragment::GetFromData(const int & A, const int & Z) {
    std::ifstream inFile("./Configs/NuclearData/nubtab12.asc");
    std::string line;
    int charge, mass;

    // read the data
    if (inFile.is_open()) {
        while (std::getline(inFile,line)) {

            charge = std::atoi(line.substr(4,4).data());
            mass   = std::atoi(line.substr(0,3).data());
            if (mass == A && charge == Z*10){
                Extract(line);
                inFile.close();
                return;
            }
        }
        throw std::runtime_error( std::string("Unable to find A = ") + A + " and Z = " + Z + "\n");
    }
    else
        throw std::runtime_error( "Unable to open file nuclear data base file\n");
}

void ReactionFragment::Extract(std::string const & line) {

    Name = line.substr(11,7);

    A = std::atoi(line.substr(0,3).data());
    Z = std::atoi(line.substr(4,4).data())/10;
    MassExcess = std::atof(line.substr(18,10).data())*UNITS::keV;

    //Lifetime
    std::string s_lt_units = line.substr(69,3);
    std::string s_LifeTime = line.substr(57,12);

    replace (s_LifeTime.begin(), s_LifeTime.end(), '*' , ' ');
    replace (s_LifeTime.begin(), s_LifeTime.end(), '<' , ' ');
    replace (s_LifeTime.begin(), s_LifeTime.end(), '>' , ' ');
    replace (s_LifeTime.begin(), s_LifeTime.end(), '&' , ' ');
    replace (s_LifeTime.begin(), s_LifeTime.end(), '#' , ' ');
    s_LifeTime.erase( std::remove_if( s_LifeTime.begin(), s_LifeTime.end(), ::isspace ), s_LifeTime.end() );
    s_lt_units.erase( std::remove_if( s_lt_units.begin(), s_lt_units.end(), ::isspace ), s_lt_units.end() );

    if(s_LifeTime=="stbl")
        LifeTime=-1;
    else
        LifeTime=std::atof(s_LifeTime.data());

    if(s_lt_units=="ms")
        LifeTime*=UNITS::ms;
    else if(s_lt_units=="us")
        LifeTime*=UNITS::us;
    else if(s_lt_units=="ns")
        LifeTime*=UNITS::ns;
    else if(s_lt_units=="ps")
        LifeTime*=UNITS::ps;
    else if(s_lt_units=="fs")
        LifeTime*=UNITS::fs;
    else if(s_lt_units=="as")
        LifeTime*=UNITS::as;
    else if(s_lt_units=="zs")
        LifeTime*=UNITS::zs;
    else if(s_lt_units=="ys")
        LifeTime*=UNITS::ys;
    else if(s_lt_units=="m")
        LifeTime*=UNITS::minutes;
    else if(s_lt_units=="h")
        LifeTime*=UNITS::hours;
    else if(s_lt_units=="d")
        LifeTime*=UNITS::days;
    else if(s_lt_units=="y")
        LifeTime*=UNITS::years;
    else if(s_lt_units=="ky")
        LifeTime*=UNITS::years*1E3;
    else if(s_lt_units=="My")
        LifeTime*=UNITS::years*1E6;
    else if(s_lt_units=="Gy")
        LifeTime*=UNITS::years*1E9;
    else if(s_lt_units=="Py")
        LifeTime*=UNITS::years*1E12;

    //Spin and parity
    std::string spinparity = line.substr(79,14);
    size_t found_p = spinparity.find("+");
    size_t found_m = spinparity.find("-");

    if (found_p != std::string::npos) Parity = "+";
    if (found_m != std::string::npos) Parity = "-";

    if (spinparity.find("0")    != std::string::npos) Spin = 0   ;
    if (spinparity.find("1")    != std::string::npos) Spin = 1   ;
    if (spinparity.find("2")    != std::string::npos) Spin = 2   ;
    if (spinparity.find("3")    != std::string::npos) Spin = 3   ;
    if (spinparity.find("4")    != std::string::npos) Spin = 4   ;
    if (spinparity.find("5")    != std::string::npos) Spin = 5   ;
    if (spinparity.find("6")    != std::string::npos) Spin = 6   ;
    if (spinparity.find("7")    != std::string::npos) Spin = 7   ;
    if (spinparity.find("8")    != std::string::npos) Spin = 8   ;
    if (spinparity.find("9")    != std::string::npos) Spin = 9   ;
    if (spinparity.find("10")   != std::string::npos) Spin = 10  ;
    if (spinparity.find("1/2")  != std::string::npos) Spin = 0.5 ;
    if (spinparity.find("3/2")  != std::string::npos) Spin = 1.5 ;
    if (spinparity.find("5/2")  != std::string::npos) Spin = 2.5 ;
    if (spinparity.find("7/2")  != std::string::npos) Spin = 3.5 ;
    if (spinparity.find("9/2")  != std::string::npos) Spin = 4.5 ;
    if (spinparity.find("11/2") != std::string::npos) Spin = 5.5 ;
    if (spinparity.find("13/2") != std::string::npos) Spin = 6.5 ;
    if (spinparity.find("15/2") != std::string::npos) Spin = 7.5 ;
    if (spinparity.find("17/2") != std::string::npos) Spin = 8.5 ;
    if (spinparity.find("19/2") != std::string::npos) Spin = 9.5 ;
    if (spinparity.find("21/2") != std::string::npos) Spin = 10.5;
}

void ReactionFragment::UpdateMass() {
    M = A * UNITS::PHYSICS::amu_c2 + MassExcess - (Z-Q) * UNITS::PHYSICS::emass_c2 + Ex;
    M2 = M * M;
}

void ReactionFragment::Set_Ex(const long double &Ex) {
    this->Ex = Ex;
    UpdateMass();
    Lab.P.SetVectM(Lab.P.Vect(), M);
    BoostToCm();
}

void ReactionFragment::Set_Ek(const long double &Ek) {
    Set_E(Ek+M);
}

void ReactionFragment::Set_Ek_cm(const long double &Ek) {
    Set_E_cm(Ek+M);
    BoostToLab();
}

void ReactionFragment::Set_E(const long double &E) {
    TVector3 tmp_vec = Lab.P.Vect();
    if (tmp_vec.Mag()==0){
        tmp_vec = sqrtl(E * E - M2)*TVector3(0., 0., 1.);
    }else {
        tmp_vec *= sqrtl(E * E - M2) / tmp_vec.Mag();
    }
    Lab.P.SetVectM(tmp_vec, M);
    BoostToCm();
}

void ReactionFragment::Set_E_cm(const long double &E) {
    TVector3 tmp_vec = Cm.P.Vect();
    if (tmp_vec.Mag()==0){
        tmp_vec = sqrtl(E * E - M2)*TVector3(0., 0., 1.);
    }else {
        tmp_vec *= sqrtl(E * E - M2) / tmp_vec.Mag();
    }
    Cm.P.SetVectM(tmp_vec, M);
    BoostToLab();
}

void ReactionFragment::Set_Pos(const TVector3 & Pos) {
    Lab.Pos = Pos;
}

void ReactionFragment::Set_Pos_cm(const TVector3 & Pos) {
    Cm.Pos = Pos;
}

void ReactionFragment::Set_Dir(const TVector3 &Dir) {
    Lab.P.Vect().SetMagThetaPhi(Lab.P.Mag(), Dir.Theta(), Dir.Phi());
    BoostToCm();
}

void ReactionFragment::Set_Dir_cm(const TVector3 &Dir) {
    Cm.P.Vect().SetMagThetaPhi(Cm.P.Mag(), Dir.Theta(), Dir.Phi());
    BoostToLab();
}

void ReactionFragment::Set_Beta(const TVector3 &Beta) {
    Set_P(M/(sqrtl(1-Beta.Mag2()))*Beta);
}

void ReactionFragment::Set_Beta_cm(const TVector3 &Beta) {
    Set_P_cm(M/(sqrtl(1-Beta.Mag2()))*Beta);
}

void ReactionFragment::Set_P(const TLorentzVector & P) {
    if (sqrt(P*P - M2)< precision) {
        Lab.P = P;
    } else {
        throw std::runtime_error(std::string("Lorentz vector does not match, computed M : ") +
                                 +sqrt(P * P)
                                 + " instead of : "
                                 +M
                                 +"\n");
    }
}

void ReactionFragment::Set_P_cm(const TLorentzVector & P) {
    if (precision > abs(P*P - M2)) {
        Cm.P = P;
    } else
        throw std::runtime_error("Lorentz vector in Reaction Fragment does not match with Fragment\n");
}

void ReactionFragment::Set_P(const TVector3 & P) {
    Lab.P.SetVectM(P, M);
    BoostToCm();
}

void ReactionFragment::Set_P_cm(const TVector3 & P) {
    Cm.P.SetVectM(P, M);
    BoostToLab();
}

void ReactionFragment::Set_P(const long double & P) {
    Lab.P.SetVectM(P/Lab.P.Vect().Mag() * Lab.P.Vect(), M);
    BoostToCm();
}

void ReactionFragment::Set_P_cm(const long double & P) {
    Cm.P.SetVectM(P/Cm.P.Vect().Mag() * Cm.P.Vect(), M);
    BoostToLab();
}

void ReactionFragment::Set_Theta_Phi(const long double & Theta, const long double & Phi) {
    TVector3 tmp_vec = Lab.P.Vect();
    tmp_vec.SetMagThetaPhi(Lab.P.Vect().Mag(), Theta, Phi);
    Lab.P.SetVectM(tmp_vec, M);
    BoostToCm();
}

void ReactionFragment::Set_Theta_Phi_cm(const long double & Theta, const long double & Phi) {
    Cm.P.Vect().SetMagThetaPhi(Cm.P.Vect().Mag(), Theta, Phi);
    BoostToLab();
}

void ReactionFragment::Set_Theta(const long double & Theta) {
    TVector3 tmp_vec = Lab.P.Vect();
    tmp_vec.SetTheta(Theta);
    Lab.P.SetVectM(tmp_vec, M);
    BoostToCm();
}

void ReactionFragment::Set_Theta_cm(const long double & Theta) {
    if (Cm.P.Vect().Mag2() == 0){
        TVector3 tmp_vec(0., 0., 1.);
        tmp_vec.SetTheta(Theta);
        Cm.P.SetVectM(tmp_vec, M);
    }else {
        TVector3 tmp_vec = Cm.P.Vect();
        tmp_vec.SetTheta(Theta);
        Cm.P.SetVectM(tmp_vec, M);
    }
    BoostToLab();
}

void ReactionFragment::Set_E_Theta(const long double & E, const long double & Theta){
    Set_E(E);
    Set_Theta(Theta);
}

void ReactionFragment::Set_E_Theta_cm(const long double & E, const long double & Theta){
    Set_Theta_cm(Theta);
    Set_E_cm(E);
}

void ReactionFragment::Set_E_Theta_Phi(const long double & E, const long double & Theta, const long double & Phi){
    Set_E(E);
    Set_Theta_Phi(Theta, Phi);
}

void ReactionFragment::Set_E_Theta_Phi_cm(const long double & E, const long double & Theta, const long double & Phi){
    Set_Theta_Phi_cm(Theta, Phi);
    Set_E_cm(E);
}

void ReactionFragment::Set_betacm(const TVector3 & beta) {
    betacm.first = true;
    betacm.second = beta;
}

void ReactionFragment::Set_Fixed(const bool & Fixed) {
    this->Fixed = Fixed;
}

void ReactionFragment::Set_ExFixed(const bool & ExFixed) {
    this->ExFixed = ExFixed;
}

bool ReactionFragment::Check_Consistent() {
    return      (precision > abs(Cm.P * Cm.P - M2))
            &&  (precision > abs(Lab.P * Lab.P - M2));
}

void ReactionFragment::BoostToCm() {
    if (betacm.first) {
        Cm.P = Lab.P;
        Cm.P.Boost(-betacm.second);
    }
}

void ReactionFragment::BoostToLab() {
    if (betacm.first) {
        Lab.P = Cm.P;
        Lab.P.Boost(betacm.second);
    }
}

long double ReactionFragment::Get_BindingEnergy() const {
    return Z*UNITS::PHYSICS::pmass_c2 + (A-Z)*UNITS::PHYSICS::nmass_c2 + (Z-Q)*UNITS::PHYSICS::emass_c2 - A * UNITS::PHYSICS::amu_c2 -  MassExcess;
}
