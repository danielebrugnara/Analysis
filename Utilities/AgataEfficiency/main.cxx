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
    for (int ii=0; ii<argc-1;++ii){
        std::string val(argv[ii+1]);
        if (val.find("--run-selector") != std::string::npos)
            run_selector = true;
    }

    for (int ii=0; ii<argc-1;++ii){
        if (run_selector) {
            ++ii;
            std::cout << "Running selector for : " << argv[ii + 1] << " \n";
            RunSelector run(argv[ii+1]);
            std::cout << "Running analysis\n";
            SpectrumAnalyzer analyze(run.GetFileName());
        }else{
            std::cout << "Running analysis\n";
            SpectrumAnalyzer analyze(argv[ii+1]);
        }
    }
    std::cout << "Finished the analysis\n";
    return 0;
}
