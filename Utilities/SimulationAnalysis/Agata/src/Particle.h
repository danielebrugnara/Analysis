#ifndef __Particle__
#define __Particle__

#include "Hit.h"
#include "TObject.h"

class Particle: public TObject{
    public:
    Particle(){};
    Particle(int type, double energy, double x, double y, double z, int nr, double beta=0, double mass=0){
        this->type=type;
        this->energy=energy;
        this->beta=beta;
        this->mass=mass;
        this->x=x;
        this->y=y;
        this->z=z;
        this->nr=nr;
    };
 //   std::vector <Hit> hits;
    void SetMass(double mass){this->mass = mass;};
    void SetEnergy(double energy){this->energy = energy;};
    private:
    int type;
//    int hits_number;
    double energy;
    double beta;
    double mass;
    double x;
    double y;
    double z;
    int nr;
    
    ClassDef(Particle, 3);
}; 

#endif
