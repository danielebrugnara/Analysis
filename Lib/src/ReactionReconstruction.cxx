#include <ReactionReconstruction.h>

ReactionReconstruction::ReactionReconstruction( ReactionFragment const & p1,
                                                ReactionFragment const & p2,
                                                ReactionFragment const & p3,
                                                ReactionFragment const & p4):
        p1(p1),
        p2(p2),
        p3(p3),
        p4(p4){}

std::string ReactionReconstruction::GetName() {
    return  p1.Get_Name()+
            "@"+
            std::to_string(p1.Get_Ek()/UNITS::MeV)+
            "("+
            p2.Get_Name()+
            ","+
            p3.Get_Name()+
            ")";
}

ReactionReconstruction2body::ReactionReconstruction2body(ReactionInput2body const & data):
        ReactionReconstruction(data[0], data[1], data[2], data[3])
{

}

ReactionReconstruction3body::ReactionReconstruction3body(ReactionInput2body const & data,
                                                         ReactionFragment::FragmentSettings const & p5):
        ReactionReconstruction2body(data),
        p5(p5){

}
