#include "ReactionReconstruction.h"
#include "ReactionFragment.h"
#include "Units.h"

#include <iostream>

int main(){
    std::cout << "Starting tester!\n";
    ReactionReconstruction2body::ReactionInput2body data{
        ReactionFragment::FragmentSettings(46, 18, 17, 10*46., 0 ),
        ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
        ReactionFragment::FragmentSettings(47, 19,17, 10*47, 0),
        ReactionFragment::FragmentSettings(2, 1, 1, 10, 0)
    };
    ReactionReconstruction2body reaction(data);
    reaction.SetBeamEnergy(600 * UNITS::MeV);
    reaction.Set_Ek(488.293*UNITS::MeV);
    const ReactionFragment & ref = reaction.GetReactionFragment(3);
    double value = ref.Get_P().Vect().Theta();
    return 0;
}
