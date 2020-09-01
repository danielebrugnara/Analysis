#ifndef EFFANALYSIS_FITTER_H
#define EFFANALYSIS_FITTER_H

#include "TH1.h"
#include "TF1.h"

#include <TFitResult.h>

class Fitter {
public:
    Fitter(TH1D spec, std::vector<std::pair<double,double>> energies);
    ~Fitter();
    struct FitRes{
        double energy;
        std::pair<double, double> integral;
        std::pair<double, double> sigma;
        std::pair<double, double> sigma_gauss;
        std::pair<double, double> tail;
        FitRes(const double& energy,
               const double& integral,
               const double& integral_err,
               const double& sigma,
               const double& sigma_err,
               const double& sigma_gauss,
               const double& sigma_gauss_err,
               const double& tail,
               const double& tail_err):
                energy(energy),
                integral(std::make_pair(integral, integral_err)),
                sigma(std::make_pair(sigma, sigma_err)),
                sigma_gauss(std::make_pair(sigma_gauss, sigma_gauss_err)),
                tail(std::make_pair(tail, tail_err)){}
    };
    std::vector<FitRes> Fit();
    TH1D& GetSpecRef(){return (TH1D &)spec;}
    TF1& GetFitRef(){return (TF1 &)fitfunc;}
private:
    TF1 fitfunc;
    TH1D spec;
    std::vector<std::pair<double, double>> energies;
};


#endif //EFFANALYSIS_FITTER_H
