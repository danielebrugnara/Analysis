#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <TSpline.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class Interpolation {
   public:
    Interpolation(std::string);
    double Evaluate(double);
    TSpline3 *GetSpline();
    double GetXPoint(int);
    double GetYPoint(int);
    ~Interpolation();

   private:
    bool ReadFile(std::string);
    std::vector<double> fValues_x;
    std::vector<double> fValues_y;
    TSpline3 *spline;
    //ClassDef(Interpolation, 1);
};

#endif