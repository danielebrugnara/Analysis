#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TSpline.h>
#include <TFile.h>
#include <TKey.h>

class Interpolation
{
public:
    explicit Interpolation(const std::string&); //Reads from mathematica txt file
    explicit Interpolation(TFile *);     //Opens existing Spline in ROOT file
    void SetUnits(const double& x_units, const double& y_units);
    inline double Evaluate(const double &x) { return spline->Eval(x*x_unit)*y_unit; };
    inline void Draw(const std::string& opt) { spline->Draw(opt.c_str()); };
    TSpline *GetSpline();
    ~Interpolation();

private:
    bool ReadFile(const std::string&);
    TSpline *spline;
    double x_unit{1};
    double y_unit{1};
    //ClassDef(Interpolation, 1);
};