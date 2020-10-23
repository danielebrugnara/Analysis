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
    reaction.ChooseFixed(3);
    reaction.ChooseExFixed(3);
    //reaction.GetReactionFragment(4).Set_Ex(0*UNITS::MeV);
    //reaction.GetReactionFragment(3).Set_Ex(25*UNITS::MeV);
    //reaction.Set_Ek((long double)(569.5*UNITS::MeV));
    //reaction.Set_Theta(0.15, true);//Works
    //std::cout << "result : " << reaction.GetFixedFragment().Get_Ek() << std::endl;//Works
    std::cout << "result : " << reaction.GetFixedFragment().Get_P().Theta() << std::endl;//Works
    //std::cout << "result : " << reaction.Get_ThetaMax()*180./3.1415 << std::endl;//works
    //reaction.Set_Theta(0.15, true);//Works
//    reaction.Set_Ek(600);
//    long double th = reaction.GetReactionFragment(3).Get_P().Vect().Theta();
//    long double a =  reaction.Set_Ek_Theta(600, th, false, true);
    //long double b =  reaction.Set_Ek_Theta(575.8 , 0.17554890526644,true, true);

    //reaction.GetReactionFragment(3).Set_Ex(10*UNITS::MeV); //Understand why!!!
//    std::cout << "result check " << a <<std::endl;
    //std::cout << "result final " << b <<std::endl;

    return 0;
}
