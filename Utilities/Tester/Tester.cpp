#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

#include <iostream>

int main(){
    std::cout << "Starting tester!\n";
    ReactionReconstruction2body::ReactionInput2body data{
        ReactionFragment::FragmentSettings(48, 20, 0, 6711.82*UNITS::MeV, 0 ),
        ReactionFragment::FragmentSettings(9, 4, 0, 0, 0*UNITS::MeV),
        ReactionFragment::FragmentSettings(48, 20,0, 0, 0*UNITS::MeV),
        ReactionFragment::FragmentSettings(9, 4, 0, 10, 0*UNITS::MeV)
    };
    ReactionReconstruction2body reaction(data);
    reaction.SetBeamEnergy(800 * UNITS::MeV);
    //reaction.Set_Ek(6232.4*UNITS::MeV);
    reaction.ChooseFixed(4);
    reaction.ChooseExFixed(3);


    reaction.Set_Theta_cm(80*UNITS::deg);
    std::cout << "deg cm: " << reaction.GetReactionFragment(4).Get_P_cm().Theta() << std::endl;
    std::cout << "deg lab: " << reaction.GetReactionFragment(4).Get_P().Theta() << std::endl;
    std::cout << "ener cm: " << reaction.GetReactionFragment(4).Get_Ek_cm() << std::endl;
    std::cout << "ener cm1: " << reaction.GetReactionFragment(1).Get_Ek_cm() << std::endl;
    std::cout << "ener cm2: " << reaction.GetReactionFragment(2).Get_Ek_cm() << std::endl;
    std::cout << "ener cm3: " << reaction.GetReactionFragment(3).Get_Ek_cm() << std::endl;
    std::cout << "ener lab: " << reaction.GetReactionFragment(4).Get_Ek() << std::endl;
    std::cout << "ener lab1: " << reaction.GetReactionFragment(1).Get_Ek() << std::endl;
    std::cout << "ener lab2: " << reaction.GetReactionFragment(2).Get_Ek() << std::endl;
    std::cout << "ener lab3: " << reaction.GetReactionFragment(3).Get_Ek() << std::endl;
    std::cout << "ener lab4: " << reaction.GetReactionFragment(4).Get_Ek() << std::endl;
    std::cout << "ang lab3: " << reaction.GetReactionFragment(3).Get_P().Theta()/UNITS::deg << std::endl;
    std::cout << "ang lab4: " << reaction.GetReactionFragment(4).Get_P().Theta()/UNITS::deg << std::endl;
    std::cout << "ang cm3: " << reaction.GetReactionFragment(3).Get_P_cm().Theta()/UNITS::deg << std::endl;
    std::cout << "ang cm4: " << reaction.GetReactionFragment(4).Get_P_cm().Theta()/UNITS::deg << std::endl;

    return 0;
}
