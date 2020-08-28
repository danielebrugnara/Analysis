#ifndef __Hit__
#define __Hit__

#include "TObject.h"

class Hit: public TObject{
    public:
        Hit(){};
        Hit(int crystal, double Energy, double x, double y, double z, int segment, int detector){
            this->crystal=crystal;
            this->Energy=Energy;
            this->x=x;
            this->y=y;
            this->z=z;
            this->segment=segment;
            this->detector=detector;
        };
        inline int GetCrystal()const{return crystal;};
        inline double GetEnergy()const{return Energy;};
        inline double GetX()const{return x;};
        inline double GetY()const{return y;};
        inline double GetZ()const{return z;};
        inline int GetSegment()const{return segment;};
        inline int GetDetector()const{return detector;};
    private:
        int crystal;
        double Energy;
        double x;
        double y;
        double z;
        int segment;
        int detector;
    
    ClassDef(Hit, 1);
};

#endif
