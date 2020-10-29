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
    reaction.ChooseExFixed(4);



    std::cout << "res : " << reaction.Set_Ek_cm(97.90) << std::endl;

    return 0;
}
