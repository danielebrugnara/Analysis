#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

#include <iostream>

int main(){
    std::cout << "Starting tester!\n";
    ReactionReconstruction2body::ReactionInput2body data{
        ReactionFragment::FragmentSettings(48, 20, 20, 6711.82*UNITS::MeV, 0 ),
        ReactionFragment::FragmentSettings(9, 4, 0, 0, 0),
        ReactionFragment::FragmentSettings(48, 20,20, 0, 0),
        ReactionFragment::FragmentSettings(9, 4, 1, 10, 0)
    };
    ReactionReconstruction2body reaction(data);
    reaction.SetBeamEnergy(800 * UNITS::MeV);
    //reaction.Set_Ek(6232.4*UNITS::MeV);
    reaction.ChooseFixed(3);
    //double en = (reaction.GetFixedFragment().Get_M())/UNITS::PHYSICS::amu_c2 * 224 * UNITS::MeV;
    reaction.Set_Ek(380.43*UNITS::MeV);
    const ReactionFragment & ref = reaction.GetReactionFragment(3);
    double value = ref.Get_P().Vect().Theta();
    return 0;
}
