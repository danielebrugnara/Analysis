#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TSpline.h>
#include <TFile.h>
#include <TKey.h>

class Interpolation {
   public:
    Interpolation(std::string);//Reads from mathematica txt file
    Interpolation(TFile *);//Opens existing Spline in ROOT file
    inline double Evaluate(double);
    TSpline *GetSpline();
    ~Interpolation();

   private:
    bool ReadFile(std::string);
    TSpline *spline;
    //ClassDef(Interpolation, 1);
};

#endif