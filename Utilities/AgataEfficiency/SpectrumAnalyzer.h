//
// Created by daniele on 13/07/2020.
//

#ifndef EFFANALYSIS_SPECTRUMANALYZER_H
#define EFFANALYSIS_SPECTRUMANALYZER_H

#include <string>

#include "DiaGraph.h"
#include "LevelScheme.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>


class SpectrumAnalyzer {
public:
    explicit SpectrumAnalyzer(const std::string &);
private:
    TH2D gg;
    void SetupLevelSchemes();

    struct Gamma{
        double egamma;
        double br;
//        static Gamma FromLevels(const Gamma & g1, const double&e1, const Gamma& g2, const double& e2){
//            return Gamma{e1-e2, g1.br};
//        }
    };

    friend std::ostream& operator<<(std::ostream& os, const Gamma& g)
    {
        os << "("<<g.egamma<<" , "<<g.br<<")";
        return os;
    }
};


#endif //EFFANALYSIS_SPECTRUMANALYZER_H
