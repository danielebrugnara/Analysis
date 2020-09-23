#include <iostream>
#include <string>
#include <memory>

#include "RunSelector.h"
#include "Readsimu.h"

#include <TFile.h>
#include <TChain.h>
#include <TVector2.h>

int main(int argc, char* argv[]){
    std::cout << "Starting analysis now \n";
    if (argc < 2)    throw std::runtime_error("set inputs correctly\n");
    std::vector<std::string> root_files;
    bool sum_all = false;
    bool is_gun = false;
    for (int ii=1; ii<argc;++ii) {
        if (std::string(argv[ii]).compare("--sum-all") == 0) {
            std::cout << "Summing all spectra in output : spectra_simu.root\n";
            sum_all = true;
            continue;
        }
        if (std::string(argv[ii]).compare("--gun") == 0) {
            std::cout << "Is gun simulation\n";
            is_gun = true;
            continue;
        }
        std::string file_name = argv[ii];
        if (file_name.compare(file_name.size() - 5, 5, ".root")) {
            bool smearing = true;
            int it_max = -1;
            std::cout << "Converting simulation for: " << file_name << std::endl;
            int nevts = readsimu(file_name, file_name + ".root", smearing, it_max, is_gun);
            file_name = file_name + ".root";
            root_files.push_back(file_name);
            std::cout << "\nnNumber of events read by readimu: " << nevts << std::endl;
        } else {
            root_files.push_back(file_name);
        }
    }

    if (sum_all) {
        TChain chain("SimulatedAgata");
        for (const auto &it: root_files) {
            std::cout << "Running selector for chain" << std::endl;
            chain.Add(it.c_str());
        }
        RunSelector run(chain, "spectra_sum_all.root");
    }else{
        for (const auto &it: root_files) {
            std::cout << "Running selector for : "<< it << std::endl;
            std::unique_ptr<TFile> file (new TFile(it.c_str()));
            if (! file->IsOpen()) throw std::runtime_error("File not valid\n");
            auto* tree = (TTree*)file->Get("SimulatedAgata");
            if (!tree) throw std::runtime_error("Tree not valid\n");

            std::string file_out_name = it;
            file_out_name.insert(0, "spectra_");
            RunSelector run(*tree, file_out_name);
            file->Close();
        }
    }
    return 0;
}
