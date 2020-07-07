#ifndef __Hit__
#define __Hit__

class Hit{
    public:
        //Hit()=delete;
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
        inline int GetSegment()const{return segment;};
        inline double GetX()const{return x;};
        inline double GetY()const{return y;};
        inline double GetZ()const{return z;};
        inline int GetDetector()const{return detector;};
        inline double GetEnergy()const{return Energy;};
    private:
        int crystal;
        int segment;
        double x;
        double y;
        double z;
        int detector;
        double Energy;
    
 //   ClassDef(Hit, 1);
};

#endif
