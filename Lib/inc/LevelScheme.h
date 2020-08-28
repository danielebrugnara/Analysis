#ifndef EFFANALYSIS_LEVELSCHEME_H
#define EFFANALYSIS_LEVELSCHEME_H

#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_map>

#include <cstring>

#include "DiaGraph.h"
#include "Units.h"



class Level{
public:
    enum Parity{
        POSITIVE_PARITY,
        NEGATIVE_PARITY,
        UNKNOWN
    };
    std::pair<double, double> Elevel;
    std::pair<double, double> T12;
    std::vector<std::pair<double, Parity>> SpinParity;
    Level(  std::pair<double, double>  Elevel,
            std::pair<double, double>  T12,
            std::vector<std::pair<double, Parity>>  SpinParity):
                Elevel(std::move(Elevel)),
                T12(std::move(T12)),
                SpinParity(std::move(SpinParity)){};
    Level()= default;
    friend std::ostream& operator<<(std::ostream& os, const Level& l)
    {
        os << "[level : "<< l.Elevel.first/UNITS::keV <<"("<<l.Elevel.second/UNITS::keV<<") keV {";
        os <<"T12 : "<< l.T12.first/UNITS::nanosecond<<"("<<l.T12.second/UNITS::nanosecond<<")}";
        return os;
    }
};

class Gamma{
public:
    std::pair<double, double> Egamma;
    std::pair<double, double> Ef_Ei;
    std::pair<int, int> Idxf_Idxi;
    std::pair<double, double> Intensity;
    std::pair<double, double> Br;
    std::string Multipol;
    Gamma(std::pair<double, double> Egamma,
          std::pair<double, double> Ef_Ei,
          std::pair<int, int> Idxf_Idxi,
          std::pair<double, double> Intensity,
          std::pair<double, double> Br,
          std::string Multipol):
                    Egamma(std::move(Egamma)),
                    Ef_Ei(std::move(Ef_Ei)),
                    Idxf_Idxi(std::move(Idxf_Idxi)),
                    Intensity(std::move(Intensity)),
                    Br(std::move(Br)),
                    Multipol(std::move(Multipol)){};
    Gamma()=default;
    friend std::ostream& operator<<(std::ostream& os, const Gamma& g)
    {
        os << "[gamma ("<<g.Egamma.first /UNITS::keV<<"("<<g.Egamma.second/UNITS::keV<<")keV) ";
        os<<": ("<<g.Ef_Ei.second/UNITS::keV <<"-->"<<g.Ef_Ei.first/UNITS::keV<<")]";
        os << "{Int : "<<g.Intensity.first <<" (" <<g.Intensity.second<<")";
        os << " , Br : "<<g.Br.first <<" (" <<g.Br.second<<")}";
        return os;
    }

};

class LevelScheme {
public:
    explicit LevelScheme(const std::string &);
    ~LevelScheme();
    std::vector<std::pair<const Gamma*, const Gamma*>> GetGammaGamma();
private:
    void ReadFile(const std::string &);
    static std::vector<std::string> ReadCSVRow(const std::string &);
    static std::pair<double, double> GetEnergyandError(const std::string &);
    static std::pair<double, double> GetLifetimeandError(const std::string &);
    static std::vector<std::pair<double, double>> GetMultipleEnergiesandError(const std::string &);
    static std::vector<std::pair<double, double>> GetMultipleIntensitiesandError(const std::string &);
    static std::vector<std::string> GetMultipolariy(const std::string &);
    static std::vector<std::pair<double, Level::Parity>> GetSpinParity(const std::string &);
    static std::vector<std::pair<double, double>> ComputeBranching(const std::vector<std::pair<double, double>> &);
    DiaGraph<Gamma, Level>* scheme;
    std::vector<std::pair<const Gamma*, const Gamma*>> gamma_gamma;
};


#endif //EFFANALYSIS_LEVELSCHEME_H
