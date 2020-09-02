#ifndef EFFANALYSIS_FITTER_H
#define EFFANALYSIS_FITTER_H

#include "TH1.h"
#include "TF1.h"

#include <unordered_map>
#include <fstream>
#include <sstream>

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
    std::vector<FitRes> Fit(const std::string&, const int&);
    TH1D& GetSpecRef(){return (TH1D &)spec;}
    TF1& GetFitRef(){return (TF1 &)fitfunc;}
private:
    void FindParameters();
    std::vector<FitRes> FindFunction();
    TF1 fitfunc;
    TH1D spec;
    std::vector<std::pair<double, double>> energies;

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
        void SetupParameter(TF1& func, int idx, const std::string& name_passed){
            this->name = name_passed;
            func.SetParameter(idx,initial_value);
            func.SetParLimits(idx,min,max);
            func.SetParName(idx, name_passed.c_str());
            if (fixed)
                func.FixParameter(idx,initial_value);
            else
                func.ReleaseParameter(idx);
        }
        void GetParameter(TF1& func, int idx){
            fitted_value.first = func.GetParameter(idx);
            fitted_value.second = func.GetParError(idx);
        }
    };

    //All parameters together
    struct InitializationParameters{
        std::vector<InitialState> ampl;     //0+ii*npars
        std::vector<InitialState> mean;     //1+ii+npars
        std::vector<InitialState> sigma;    //2+ii*npars
        std::vector<InitialState> tau;      //3+ii*npars
        InitialState common_ampl;           //ampl.size()*npars
        InitialState common_offset;         //ampl.size()*npars+1

        double left_interval{};
        double right_interval{};
        double min_sigma{};
        double max_sigma{};
        double min_tau{};
        double max_tau{};

        int npars{4};
        explicit InitializationParameters(int n_gauss){
            ampl.resize(n_gauss);
            mean.resize(n_gauss);
            sigma.resize(n_gauss);
            tau.resize(n_gauss);
        };
        explicit InitializationParameters()= default;
        InitializationParameters(const InitializationParameters& old_obj)=default;
        void SetupParameters(TF1& func){
            for (unsigned long int ii=0; ii<ampl.size(); ++ii){
                ampl[ii].SetupParameter(func,0+ii*npars,Form("Ampl_%li", ii));
                mean[ii].SetupParameter(func,1+ii*npars,Form("Mean_%li", ii));
                sigma[ii].SetupParameter(func,2+ii*npars,Form("Sigma_%li", ii));
                tau[ii].SetupParameter(func,3+ii*npars,Form("Tau_%li", ii));
            }
            common_ampl.SetupParameter(func,(int)ampl.size()*npars,"CommonAmpl");
            common_offset.SetupParameter(func,(int)ampl.size()*npars+1,"CommonOffset");
            func.SetRange(left_interval, right_interval);
        }
        void GetParameters(TF1& func){
            for (unsigned long int ii=0; ii<ampl.size(); ++ii){
                ampl[ii].GetParameter(func,0+(int)ii*npars);
                mean[ii].GetParameter(func,1+(int)ii*npars);
                sigma[ii].GetParameter(func,2+(int)ii*npars);
                tau[ii].GetParameter(func,3+(int)ii*npars);
            }
            common_ampl.GetParameter(func, (int)ampl.size()*npars);
            common_offset.GetParameter(func, (int)ampl.size()*npars+1);
        }
    };
    InitializationParameters parameters;
public:
    void WriteParsOnFile(const std::string&, const int&);
private:
    InitializationParameters ReadParsFromFile(const std::string&, const int&);
};


#endif //EFFANALYSIS_FITTER_H
