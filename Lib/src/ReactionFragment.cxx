#include <ReactionFragment.h>

ReactionFragment::ReactionFragment(FragmentSettings const & data):
                                            A   (data.A),
                                            Z   (data.Z),
                                            Q   (data.Q),
                                            M   (0),
                                            Ek  (data.Ek),
                                            E   (0),
                                            Ex  (data.Ex),
                                            P   (TLorentzVector()),
                                            Pos (TVector3())
{
    GetFromData(A, Z);
    this->M = A*UNITS::PHYSICS::amu_c2 + MassExcess - Q*UNITS::PHYSICS::emass_c2+Ex;
    this->E = Ek + this->M;
    P.SetPxPyPzE(0., 0.,sqrt(pow(E, 2) - pow(M, 2)), E);
    P.SetVect(data.Dir);
}

void ReactionFragment::GetFromData(const int & A, const int & Z) {
    std::ifstream inFile("./Configs/NuclearData/nubtab12.asc");
    std::string line;
    int charge, mass;

    // read the data
    if (inFile.is_open()) {
        while (std::getline(inFile,line)) {

            charge = atoi(line.substr(4,4).data());
            mass   = atoi(line.substr(0,3).data());
            if (mass == A && charge == Z*10){
                Extract(line.data());
                inFile.close();
                return;
            }
        }
        throw std::runtime_error( std::string("Unable to find A = ") + A + " and Z = " + Z + "\n");
    }
    else
        throw std::runtime_error( "Unable to open file nuclear data base file\n");
}

void ReactionFragment::Extract(std::string line) {

    Name = line.substr(11,7);

    A = atoi(line.substr(0,3).data());
    Z = atoi(line.substr(4,4).data())/10;
    MassExcess = atof(line.substr(18,10).data())*UNITS::keV;

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
        LifeTime=atof(s_LifeTime.data());

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
    std::string spinparity = line.substr(79,14).data();
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
