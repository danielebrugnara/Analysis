#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TSpline.h>
#include <TGraph.h>
#include <TProfile.h>
#include <TFile.h>
#include <TKey.h>

class Interpolation
{
public:
    explicit Interpolation(const std::string& ); //Reads from mathematica txt file
    explicit Interpolation(TFile *, const bool& useHistogram=false, const bool& useGraph=false);     //Opens existing Spline in ROOT file
    void SetUnits(const double& x_units, const double& y_units);
    inline double Evaluate(const double &x) {
        if(spline!= nullptr) return spline->Eval(x*x_unit)*y_unit;
        if(hist!= nullptr) return hist->GetBinContent( hist->FindBin(x*x_unit))*y_unit;
        if(graph!= nullptr) return graph->Eval(x*x_unit)*y_unit;
        throw std::runtime_error("Interpolations are nullptr\n");
    };
    inline void Draw(const std::string& opt) { spline->Draw(opt.c_str()); };
    TSpline *GetSpline();
    ~Interpolation();

private:
    bool ReadFile(const std::string&);
    TSpline* spline{nullptr};
    TGraph* graph{nullptr};
    TProfile* hist{nullptr};
    double x_unit{1};
    double y_unit{1};
    //ClassDef(Interpolation, 1);
};