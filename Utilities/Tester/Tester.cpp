#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

#include <iostream>

int main(){
    std::cout << "Starting tester!\n";
    ReactionReconstruction2body<long  double>::ReactionInput2body data{
        ReactionFragment::FragmentSettings(46, 18, 0, 450*UNITS::MeV, 0*UNITS::MeV ),
//ReactionFragment::FragmentSettings(46, 18, 0, 373.7*UNITS::MeV, 0*UNITS::MeV ),
        ReactionFragment::FragmentSettings(3, 2, 0, 0, 0*UNITS::MeV),
        ReactionFragment::FragmentSettings(2, 1, 0, 10, 0*UNITS::MeV),
        ReactionFragment::FragmentSettings(47, 19,0, 0, 0*UNITS::MeV)
    };
    ReactionReconstruction2body<long double> reaction(data);
    reaction.ChooseFixed(3);
    reaction.ChooseExFixed(3);

    //std::cout << "Angles\n";
    //std::cout << "angle fixed cm: " << reaction.GetFixedFragment().Get_P_cm().Theta()*TMath::RadToDeg() << std::endl;
    //std::cout << "angle free cm: " << reaction.GetFreeFragment().Get_P_cm().Theta()*TMath::RadToDeg() << std::endl;
    //std::cout << "angle fixed lab: " << reaction.GetFixedFragment().Get_P().Theta()*TMath::RadToDeg() << std::endl;
    //std::cout << "angle free lab: " << reaction.GetFreeFragment().Get_P().Theta()*TMath::RadToDeg() << std::endl;

    std::cout << "Energy\n";
    reaction.GetFreeFragment().Set_Ex(0*UNITS::MeV);
    reaction.Set_Theta_cm(160*UNITS::deg);
    std::cout << "theta lab fixed: " << reaction.GetFixedFragment().Get_P().Theta()*TMath::RadToDeg()<< std::endl;

    //std::cout << "energy lab fixed cm: " << reaction.GetFixedFragment().Get_Ek_cm()/UNITS::MeV << std::endl;
    //std::cout << "energy lab ree cm: " << reaction.GetFreeFragment().Get_Ek_cm()/UNITS::MeV << std::endl;
    //std::cout << "energy lab fixed : " << reaction.GetFixedFragment().Get_Ek()/UNITS::MeV << std::endl;
    //std::cout << "energy lab free : " << reaction.GetFreeFragment().Get_Ek()/UNITS::MeV << std::endl;

    std::cout << "Energy ex=2\n";
    reaction.GetFreeFragment().Set_Ex(2*UNITS::MeV);
    reaction.Set_Theta_cm(170*UNITS::deg);
    //std::cout << "energy lab fixed cm: " << reaction.GetFixedFragment().Get_Ek_cm()/UNITS::MeV << std::endl;
    //std::cout << "energy lab ree cm: " << reaction.GetFreeFragment().Get_Ek_cm()/UNITS::MeV << std::endl;
    //std::cout << "energy lab fixed : " << reaction.GetFixedFragment().Get_Ek()/UNITS::MeV << std::endl;
    //std::cout << "energy lab free : " << reaction.GetFreeFragment().Get_Ek()/UNITS::MeV << std::endl;
    return 0;
}
