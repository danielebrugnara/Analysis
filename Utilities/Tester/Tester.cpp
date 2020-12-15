#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

#include <iostream>

int main(){
    std::cout << "Starting tester!\n";
    ReactionReconstruction2body<long  double>::ReactionInput2body data{
        ReactionFragment::FragmentSettings(46, 18, 0, 458*UNITS::MeV, 0 ),
        ReactionFragment::FragmentSettings(3, 2, 0, 0, 0*UNITS::MeV),
        ReactionFragment::FragmentSettings(47, 19,0, 0, 0*UNITS::MeV),
        ReactionFragment::FragmentSettings(2, 1, 0, 10, 0*UNITS::MeV)
    };
    ReactionReconstruction2body<long double> reaction(data);
    //reaction.SetBeamEnergy(800 * UNITS::MeV);
    //reaction.Set_Ek(6232.4*UNITS::MeV);
    reaction.ChooseFixed(4);
    reaction.ChooseExFixed(4);


    reaction.Set_Theta(130*UNITS::deg,false);

    std::cout << "res : " << reaction.GetFixedFragment().Get_Ek()/UNITS::MeV << std::endl;



    reaction.Set_Theta(133*UNITS::deg,false);
    std::cout << "res : " << reaction.GetFixedFragment().Get_Ek()/UNITS::MeV << std::endl;
    return 0;
}
