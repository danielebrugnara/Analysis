#include <iostream>

#include "TApplication.h"
#include "TCanvas.h"

#include "RunSelector.h"
#include "SpectrumAnalyzer.h"

#include "Globals.h"

int main(int argc, char* argv[]){

    TApplication theApp("app", new int(0), new char*);
    std::cout << "Starting analysis now \n";
    if (argc < 2)    throw std::runtime_error("set inputs correctly\n");
    bool run_selector = false;
    bool debug_canvas = false;
    bool use_cores = false;
    //std::map<int, bool> skip_argc;
    std::vector<std::string> files;
    for (int ii=1; ii<argc;++ii){
        std::string val(argv[ii]);
        if (val.find("-h") != std::string::npos) {
            std::cout << "--run-selector\t\t\t\truns selector starting from tree\n";
            std::cout << "--debug-canvas\t\t\t\tshows canvas with fits\n";
            std::cout << "--use-cores\t\t\t\tuses core graphs instead of addback\n";
            return 0;
        }
        if (val.find("--run-selector") != std::string::npos) {
            run_selector = true;
            //skip_argc[ii] = true;
            continue;
        }
        if (val.find("--debug-canvas") != std::string::npos) {
            debug_canvas = true;
            //skip_argc[ii] = true;
            continue;
        }
        if (val.find("--use-cores") != std::string::npos) {
            use_cores = true;
            //skip_argc[ii] = true;
            continue;
        }
        files.push_back(val);
    }

    for (const auto& file:files){
        if (run_selector) {
            std::cout << "Running selector for : " << file << " \n";
            RunSelector run(file);
            std::cout << "Running analysis\n";
            SpectrumAnalyzer analyze(run.GetFileName(), debug_canvas, true);
        }else{
            std::cout << "Running analysis, file : " << file << "\n";
            SpectrumAnalyzer analyze(file, debug_canvas, use_cores);
        }
    }
    std::cout << "Finished the analysis\n";
    return 0;
}
