#ifndef __REACTIONRECONSTRUCTION_H__
#define __REACTIONRECONSTRUCTION_H__

#include <string>

#include <TLorentzVector.h>

#include <ReactionFragment.h>

class ReactionReconstruction{
protected:
    ReactionReconstruction( ReactionFragment const &,
                            ReactionFragment const &,
                            ReactionFragment const &,
                            ReactionFragment const &);
    ~ReactionReconstruction()=default;
    ReactionFragment p1;
    ReactionFragment p2;
    ReactionFragment p3;
    ReactionFragment p4;

    virtual std::string GetName();

private:
    ReactionReconstruction()=delete;

public:
    //Nothing public
};

//p1 + p2 -> p3 + p4 ## where p2 is at rest
class ReactionReconstruction2body : private ReactionReconstruction{
private:
    ReactionReconstruction2body()=delete;

public:
    typedef std::array<ReactionFragment::FragmentSettings, 4> ReactionInput2body;
    ReactionReconstruction2body(ReactionInput2body const &);
    ~ReactionReconstruction2body()=default;
};

//p1 + p2 -> p3 + p4* -> p3 + p4 + p5
class ReactionReconstruction3body : private ReactionReconstruction2body{
private:
    ReactionReconstruction3body()=delete;

public:
    ReactionReconstruction3body(ReactionInput2body const &,
                                ReactionFragment::FragmentSettings const &);
    ~ReactionReconstruction3body()=default;

private:
    ReactionFragment p5;
};

#endif