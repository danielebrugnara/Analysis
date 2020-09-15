#include <iostream>
#include <string.h>

#include "RunSelector.h"

#include <TFile.h>
#include <TVector2.h>

int main(int argc, char* argv[]){
    std::cout << "Starting analysis now \n";
    if (argc < 2)    throw std::runtime_error("set inputs correctly\n");
    int nevts = 0;
    bool sum_all = false;
    std::string command;
    command += "hadd -f spectra_simu.root ";
    for (int ii=0; ii<argc-1;++ii){
        if (std::string(argv[ii+1]).compare("--sum-all") == 0){
            std:: cout << "Summing all spectra in output : spectra_simu.root\n";
            sum_all = true;
            continue;
        }
        std::cout << "Running selector for : "<<argv[ii+1]<<" \n";
        RunSelector run(argv[ii+1]);
        std::cout << "Partial number of events : "<< run.nevts << std::endl;
        nevts += run.nevts;
        command += run.file_out_name;
        command += " ";
    }
    std::cout << "Total number of events : " << nevts << std::endl;
    if (sum_all){
        system(command.c_str());
        TFile file("spectra_simu.root", "update");
        TVector2 start_stop(0.,nevts);
        start_stop.Write("start_stop");
        file.Close();
    }
    return 0;
}
