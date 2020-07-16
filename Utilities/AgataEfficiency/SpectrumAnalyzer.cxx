#include "SpectrumAnalyzer.h"

SpectrumAnalyzer::SpectrumAnalyzer(const std::string & file_name) {
    auto* file = new TFile(file_name.c_str());
    gg = *((TH2D*) file->Get("mgamma_gamma"));

    file->Close();
    SetupLevelSchemes();

}

void SpectrumAnalyzer::SetupLevelSchemes() {
    LevelScheme("adoptedLevels152Sm.csv");
    //152Sm
//    std::vector<double> levels_152Sm = {0, 121.8, 366.5, 1085.5, 1233.9, 1529.8};
//    std::vector<graphEdge<Gamma>> gammas_152Sm = {
//            // (x, y, w) -> edge from x to y with weight w
//            {5,1,{0, 0.5}},
//            {5,3,{0, 0.5}},
//            {4,1,{0, 1}},
//            {3,0,{0, 0.5}},
//            {3,1,{0, 0.5}},
//            {2,1,{0, 1}},
//            {1,0,{0, 1}},
//            {0,0,{0, 1}}
//    };
//
//    DiaGraph<Gamma, double> scheme_152Sm(   &gammas_152Sm[0],
//                                            &levels_152Sm[0],
//                                            Gamma::FromLevels,
//                                            gammas_152Sm.size(),
//                                            levels_152Sm.size());
//    scheme_152Sm.DumpGraph();
//    //152Gd
//    std::vector<double> levels_152Gd = {0, 344.3, 1123.2};
//    std::vector<graphEdge<double>> gammas_152Gd = {
//            // (x, y, w) -> edge from x to y with weight w
//            {2,1,0},
//            {1,0,0},
//            {0,0,0}
//    };
//    DiaGraph<double> scheme_152Gd(&gammas_152Gd[0],
//                                  &levels_152Gd[0],
//                                  [](double a, double b){return a-b;},
//                                  gammas_152Gd.size(),
//                                  levels_152Gd.size());
//    scheme_152Gd.DumpGraph();
}
