#pragma once

#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussModel.h"
#include "RooPolynomial.h"
#include "RooDecay.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"


#include <unordered_map>
#include <fstream>
#include <sstream>

#include <TFitResult.h>

class Fitter {
public:
    Fitter(const TH1D&, std::vector<std::pair<double,double>>, const bool& left_tail=true);
    Fitter()=default;
    void EnableCanvas(bool canvas_enabled=true);

public:
    void SetMeanIntervalAroundCenter(const double& min, const double&max)
                    {mean_around_center_min=min; mean_around_center_max=max;};
    void SetInitialSigma(const double& val){initial_sigma=val;};
    void SetSigmaInterval(const double& min, const double&max){sigma_min=min; sigma_max=max;};
    void SetInitialTau(const double& val){initial_tau=val;};
    void SetTauInterval(const double& min, const double&max){tau_min=min; tau_max=max;};

private:
    double mean_around_center_min, mean_around_center_max;
    double initial_sigma;
    double sigma_min, sigma_max;
    double initial_tau;
    double tau_min, tau_max;


    void SetFittingInterval(const double& min, const double& max){fitting_min=min; fitting_max=max;};
    double fitting_min, fitting_max;
public:
    ~Fitter();
    struct FitRes{
        double energy;
        std::pair<double, double> integral;
        std::pair<double, double> ampl;
        std::pair<double, double> mean_gauss;
        std::pair<double, double> sigma_gauss;
        std::pair<double, double> tail;
        FitRes(const double& energy,
               const double& integral,
               const double& integral_err,
               const double& ampl,
               const double& ampl_err,
               const double& mean_gauss,
               const double& mean_gauss_err,
               const double& sigma_gauss,
               const double& sigma_gauss_err,
               const double& tail,
               const double& tail_err):
                energy(energy),
                integral(std::make_pair(integral, integral_err)),
                ampl(std::make_pair(ampl, ampl_err)),
                mean_gauss(std::make_pair(mean_gauss, mean_gauss_err)),
                sigma_gauss(std::make_pair(sigma_gauss, sigma_gauss_err)),
                tail(std::make_pair(tail, tail_err)){}
    };
    std::vector<FitRes> Fit();
    TH1D& GetSpecRef(){return (TH1D &)spec;}
    TF1& GetFitRef(){return (TF1 &)fitfunc;}
private:
    void FindParameters();
    TF1 fitfunc;
    TH1D spec;
    bool left_tail;
    std::vector<std::pair<double, double>> energies;
    bool canvas_enabled;

    //Single parameter values
    struct InitialState{
        double initial_value;
        std::pair<double, double> fitted_value;
        double min;
        double max;
        bool fixed;
        std::string name;
        InitialState(double initial_value, double min, double max, bool fixed):
                initial_value(initial_value),
                fitted_value({0,0}),
                min(min),
                max(max),
                fixed(fixed){}
        InitialState():
                initial_value(0),
                fitted_value({0,0}),
                min(0),
                max(0),
                fixed(false){}
        void GetParameter(const RooRealVar* const par) {
            fitted_value.first = par->getVal();
            fitted_value.second = par->getError();
        }

    };

    struct RooFitParams{
        RooRealVar* fracsignal;
        RooRealVar* mean;
        RooRealVar* sigma;
        RooRealVar* tau;
        RooFitParams(	RooRealVar* fracsignal,
                         RooRealVar* mean,
                         RooRealVar* sigma,
                         RooRealVar* tau):
                fracsignal(fracsignal),
                mean(mean),
                sigma(sigma),
                tau(tau){};
        RooFitParams():
                fracsignal(nullptr),
                mean(nullptr),
                sigma(nullptr),
                tau(nullptr){};
    };

    //All parameters together
    struct InitializationParameters{
        std::vector<InitialState> ampl;
        std::vector<InitialState> mean;
        std::vector<InitialState> sigma;
        std::vector<InitialState> tau;
        InitialState common_ampl;
        InitialState common_offset;

        double left_interval{};
        double right_interval{};
        double min_sigma{};
        double max_sigma{};
        double min_tau{};
        double max_tau{};

        explicit InitializationParameters(int n_gauss){
            ampl.resize(n_gauss);
            mean.resize(n_gauss);
            sigma.resize(n_gauss);
            tau.resize(n_gauss);
        };
        explicit InitializationParameters()= default;
        InitializationParameters(const InitializationParameters& old_obj)=default;
        void GetParameters(std::vector<RooFitParams>& params){
            for (unsigned long int ii=0; ii<ampl.size(); ++ii){
                ampl[ii].GetParameter(params[ii].fracsignal);
                mean[ii].GetParameter(params[ii].mean);
                sigma[ii].GetParameter(params[ii].sigma);
                tau[ii].GetParameter(params[ii].tau);
            }
        }
    };
    InitializationParameters parameters;
    bool mean_fixed;
    bool relative_amplitude_fixed;
public:
    void SetMeanFixed(const bool& mean_fixed = true)
            { this->mean_fixed = mean_fixed;};
    void SetRelativeAmplitudeFixed(const bool& relative_amplitude_fixed = true)
            { this->relative_amplitude_fixed = relative_amplitude_fixed;};
    void WriteParsOnFile(const std::string&, const int&);
    void ReadParsFromFile(const std::string&, const int&);
private:
};