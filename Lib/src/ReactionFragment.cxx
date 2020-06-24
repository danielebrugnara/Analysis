#include <ReactionFragment.h>

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
        Invariant   (-1),
        Lab         (data.Ek),
        Cm          (0),
        Fixed       (true),
        ExFixed     (true),
        precision   (1E-12)
{
    GetFromData(A, Z);
    Lab.modified = true;
    ComputeChanges();
    //P.SetPxPyPzE(0., 0.,sqrt(pow(E, 2) - M2), E);
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

ReactionFragment::Properties & ReactionFragment::GetModified() {
    if (Lab.modified && !Cm.modified)
        return Lab;
    if (!Lab.modified && Cm.modified)
        return Cm;
    else
        throw std::runtime_error("Reaction fragment not setup correctly\n");
}

ReactionFragment::Properties & ReactionFragment::GetNotModified() {
    if (!Lab.modified && Cm.modified)
        return Lab;
    if (Lab.modified && !Cm.modified)
        return Cm;
    else
        throw std::runtime_error("Reaction fragment not setup correctly\n");
}

void ReactionFragment::ComputeChanges() {
    M = A*UNITS::PHYSICS::amu_c2 + MassExcess - Q*UNITS::PHYSICS::emass_c2+Ex;
    M2 = M*M;
    auto & changed      = GetModified();
    auto & notchanged   = GetNotModified();

    changed.E = changed.Ek + M;
    TLorentzVector tmp_vec = TLorentzVector(0., 0., 0., M);
    if ((Lab.P.Vect().Mag() - 0.) < precision) {
        changed.P -= tmp_vec;
        changed.P.Boost(TVector3(0., 0., sqrt(changed.E * changed.E - M2)));
    } else {
        tmp_vec.Boost((sqrt(changed.E * changed.E - M2) / changed.P.Vect().Mag()) * changed.P.Vect());
        changed.P = tmp_vec;
    }
    if (changed.P_tot.first) {
        notchanged.P.Boost(-changed.P_tot.second.BoostVector());
        notchanged.E = notchanged.P.E();
        notchanged.Ek = notchanged.E - M;
    }
}

void ReactionFragment::Set_Ex(const double &Ex) {
    this->Ex = Ex;
    ComputeChanges();
}

void ReactionFragment::Set_Ek(const double &Ek) {
    //Renormalization of beta computed with E = gamma m c^2
    Cm.modified     = false;
    Lab.modified    = true;
    Lab.Ek = Ek;
    ComputeChanges();
}

void ReactionFragment::Set_Ek_cm(const double &Ek) {
    //Renormalization of beta computed with E = gamma m c^2
    Cm.modified     = true;
    Lab.modified    = false;
    Cm.Ek = Ek;
    ComputeChanges();
}

void ReactionFragment::Set_E(const double &E) {
    Cm.modified     = false;
    Lab.modified    = true;
    Lab.Ek = E - M;
    ComputeChanges();
}

void ReactionFragment::Set_E_cm(const double &E) {
    Cm.modified     = true;
    Lab.modified    = false;
    Cm.Ek = E - M;
    ComputeChanges();
}

void ReactionFragment::Set_Pos(const TVector3 & Pos) {
    Lab.Pos = Pos;
    //No need to compute changes
}

void ReactionFragment::Set_Pos_cm(const TVector3 & Pos) {
    Cm.Pos = Pos;
    //No need to compute changes
}

void ReactionFragment::Set_Dir(const TVector3 &Dir) {
    Cm.modified     = false;
    Lab.modified    = true;
    Lab.P.SetVect(Dir); //The correct vector will be computed anyway
    ComputeChanges();
}

void ReactionFragment::Set_Dir_cm(const TVector3 &Dir) {
    Cm.modified     = true;
    Lab.modified    = false;
    Cm.P.SetVect(Dir); //The correct vector will be computed anyway
    ComputeChanges();
}

void ReactionFragment::Set_Beta(const TVector3 &Beta) {
    Cm.modified     = false;
    Lab.modified    = true;
    Lab.P.SetVect(Beta); //The correct vector will be set based on the energy
    Lab.Ek = M / sqrt(1- Beta.Mag2()) - M;
    ComputeChanges();
}

void ReactionFragment::Set_Beta_cm(const TVector3 &Beta) {
    Cm.modified     = true;
    Lab.modified    = false;
    Cm.P.SetVect(Beta); //The correct vector will be set based on the energy
    Cm.Ek = M / sqrt(1- Beta.Mag2()) - M;
    ComputeChanges();
}

void ReactionFragment::Set_P(const TLorentzVector & P) {
    //Mass is set, so boost vector should be sufficient
    Cm.modified     = false;
    Lab.modified    = true;
    if (precision > abs(P*P - M2)) {
        Lab.P = P;
        Lab.Ek = P.E() - M;
    } else
        throw std::runtime_error(std::string("Lorentz vector does not match, computed M : ")+
                                 +(P*P)
                                 +"\n");
    ComputeChanges();
}

void ReactionFragment::Set_P_cm(const TLorentzVector & P) {
    //Mass is set, so boost vector should be sufficient
    Cm.modified     = true;
    Lab.modified    = false;
    if (precision > abs(P*P - M2)) {
        Cm.P = P;
        Cm.Ek = P.E() - M;
    } else
        throw std::runtime_error("Lorentz vector in Reaction Fragment does not match with Fragment\n");
    ComputeChanges();
}

void ReactionFragment::Set_P(const TVector3 & P) {
    Cm.modified     = false;
    Lab.modified    = true;
    Lab.P.SetVect(P);
    Lab.Ek = sqrt(P.Mag2()+M2) - M;
    ComputeChanges();
}

void ReactionFragment::Set_P_cm(const TVector3 & P) {
    Cm.modified     = true;
    Lab.modified    = false;
    Cm.P.SetVect(P);
    Cm.Ek = sqrt(P.Mag2()+M2) - M;
    ComputeChanges();
}

void ReactionFragment::Set_P(const double & P) {
    Cm.modified     = false;
    Lab.modified    = true;
    Lab.Ek = sqrt(P*P+M2) - M;
    ComputeChanges();
}

void ReactionFragment::Set_P_cm(const double & P) {
    Cm.modified     = true;
    Lab.modified    = false;
    Cm.Ek = sqrt(P*P+M2) - M;
    ComputeChanges();
}

void ReactionFragment::Set_Theta_Phi(const double & Theta, const double & Phi) {
    Cm.modified     = false;
    Lab.modified    = true;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Phi);
    Lab.P.SetVect(tmp_vec);
    ComputeChanges();
}

void ReactionFragment::Set_Theta_Phi_cm(const double & Theta, const double & Phi) {
    Cm.modified     = true;
    Lab.modified    = false;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Phi);
    Cm.P.SetVect(tmp_vec);
    ComputeChanges();
}

void ReactionFragment::Set_Theta(const double & Theta) {
    Cm.modified     = false;
    Lab.modified    = true;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Lab.P.Vect().Phi());
    Lab.P.SetVect(tmp_vec);
    ComputeChanges();
}

void ReactionFragment::Set_Theta_cm(const double & Theta) {
    Cm.modified     = true;
    Lab.modified    = false;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Cm.P.Vect().Phi());
    Cm.P.SetVect(tmp_vec);
    ComputeChanges();
}

void ReactionFragment::Set_E_Theta(const double & E, const double & Theta){
    Cm.modified     = false;
    Lab.modified    = true;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Lab.P.Vect().Phi());
    Lab.P.SetVect(tmp_vec);
    Lab.Ek = E - M;
    ComputeChanges();
}

void ReactionFragment::Set_E_Theta_cm(const double & E, const double & Theta){
    Cm.modified     = true;
    Lab.modified    = false;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Cm.P.Vect().Phi());
    Cm.P.SetVect(tmp_vec);
    Cm.Ek = E - M;
    ComputeChanges();
}

void ReactionFragment::Set_E_Theta_Phi(const double & E, const double & Theta, const double & Phi){
    Cm.modified     = false;
    Lab.modified    = true;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Phi);
    Lab.P.SetVect(tmp_vec);
    Lab.Ek = E - M;
    ComputeChanges();
}

void ReactionFragment::Set_E_Theta_Phi_cm(const double & E, const double & Theta, const double & Phi){
    Cm.modified     = true;
    Lab.modified    = false;
    TVector3 tmp_vec;
    tmp_vec.SetMagThetaPhi(1., Theta, Phi);
    Cm.P.SetVect(tmp_vec);
    Cm.Ek = E - M;
    ComputeChanges();
}

void ReactionFragment::Set_P_tot(const TLorentzVector & P_tot) {
    Cm.modified     = false;
    Lab.modified    = true;
    Lab.P_tot.first = true;
    Lab.P_tot.second = P_tot;
    ComputeChanges();
}

void ReactionFragment::Set_P_tot_cm(const TLorentzVector & P_tot) {
    Cm.modified     = true;
    Lab.modified    = false;
    Cm.P_tot.first = true;
    Cm.P_tot.second = P_tot;
    ComputeChanges();
}

double ReactionFragment::Compute_Ex(const TLorentzVector & P) {
    return P*P-M;
}

void ReactionFragment::Set_Fixed(const bool & Fixed) {
    this->Fixed = Fixed;
}

void ReactionFragment::Set_ExFixed(const bool & ExFixed) {
    this->ExFixed = ExFixed;
}

bool ReactionFragment::Check_Consistent() {
    return precision > abs(GetModified().P * GetModified().P - M2);
}
