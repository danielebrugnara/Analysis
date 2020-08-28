#ifndef __Particle__
#define __Particle__

#include "Hit.h"
#include "TObject.h"

class Particle: public TObject{
    public:
    Particle(){};
    Particle(int type, double Energy, double x, double y, double z, int nr){
        this->type=type;
        this->Energy=Energy;
        this->x=x;
        this->y=y;
        this->z=z;
        this->nr=nr;
    };
 //   std::vector <Hit> hits;
    private:
    int type;
//    int hits_number;
    double Energy;
    double x;
    double y;
    double z;
    int nr;
    
    ClassDef(Particle, 1);
}; 

#endif
