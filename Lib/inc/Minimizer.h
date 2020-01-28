#ifndef __MINIMIZER_H__
#define __MINIMIZER_H__

#include <array>

class Minimizer{
    public:
        Minimizer(double, bool, double,  double, double, double, int, double);
        ~Minimizer();

    private:
        std::array <double, 2> Y;
        std::array <double, 2> X;
        std::array < long double, 2> Derivative;
        std::array <long double, 2> Step;
        std::array <long double, 2> Rate;
        int     NSteps;
        bool    Verbose;
        double  StartingRate;
        double  Coefficient;
        double  Threshold;
        double  Quenching;
        double  MaxSteps;

};

#endif