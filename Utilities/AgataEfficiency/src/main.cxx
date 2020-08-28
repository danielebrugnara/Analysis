#include <iostream>

#include "TApplication.h"
#include "TCanvas.h"

#include "RunSelector.h"
#include "SpectrumAnalyzer.h"

#include "Globals.h"

int main(int argc, char* argv[]){

    TApplication theApp("app", new int, new char*);
    std::cout << "Starting analysis now \n";
    if (argc < 2)    throw std::runtime_error("set inputs correctly\n");
    bool run_selector = false;
    bool debug_canvas = false;
    std::map<int, bool> skip_argc;
    for (int ii=0; ii<argc-1;++ii){
        std::string val(argv[ii+1]);
        if (val.find("--run-selector") != std::string::npos) {
            run_selector = true;
            skip_argc[ii] = true;
            continue;
        }
        if (val.find("--debug-canvas") != std::string::npos) {
            debug_canvas = true;
            skip_argc[ii] = true;
            continue;
        }
        skip_argc[ii] = false;
    }

    for (int ii=0; ii<argc-1;++ii){
        if (skip_argc[ii]) continue;
        if (run_selector) {
            std::cout << "Running selector for : " << argv[ii + 1] << " \n";
            RunSelector run(argv[ii+1]);
            std::cout << "Running analysis\n";
            SpectrumAnalyzer analyze(run.GetFileName(), debug_canvas);
        }else{
            std::cout << "Running analysis, file : " << argv[ii+1] << "\n";
            SpectrumAnalyzer analyze(argv[ii+1], debug_canvas);
        }
    }
    std::cout << "Finished the analysis\n";
    return 0;
}
