#include <iostream>
#include <string>
#include <memory>

#include "RunSelector.h"
#include "Readsimu.h"

#include <TFile.h>
#include <TChain.h>
#include <TVector2.h>
#include <TApplication.h>

int main(int argc, char* argv[]){
    TApplication theApp("app", new int(0), new char*);
    std::cout << "Starting analysis now \n";
    if (argc < 2)    throw std::runtime_error("set inputs correctly\n");
    std::vector<std::string> root_files;
    bool sum_all = false;
    bool is_gun = false;
    bool infinite_resolution = false;
    std::vector<double> fit_energies;
    for (unsigned int ii=1; ii<(unsigned int)argc;++ii) {
        if (std::string(argv[ii]) == "--sum-all") {
            std::cout << "Summing all spectra in output : spectra_simu.root\n";
            sum_all = true;
            continue;
        }
        if (std::string(argv[ii]) == "--gun") {
            std::cout << "Is gun simulation\n";
            is_gun = true;
            continue;
        }
        if (std::string(argv[ii]) == "--infinite-resolution") {
            std::cout << "Has infinite resolution\n";
            infinite_resolution = true;
            continue;
        }
        if (std::string(argv[ii]) == "--fit-efficiency"){
            try {
                fit_energies.resize(std::stoi(argv[++ii]));
                std::cout << "Fitting n=" << fit_energies.size() << " energies\n";
                for(double & fit_energy : fit_energies){
                    fit_energy = std::stod(argv[++ii]);
                    std::cout << "Fitting " << fit_energy << std::endl;
                }
            } catch (...) {
                throw std::runtime_error("Error parsing fit-efficiency");
            }
            continue;
        }
        std::string file_name = argv[ii];
        if (file_name.compare(file_name.size() - 5, 5, ".root")) {
            bool smearing = !infinite_resolution;
            int it_max = -1;
            std::cout << "Converting simulation for: " << file_name << std::endl;
            int nevts = readsimu(file_name, file_name + ".root", smearing, it_max, is_gun);
            file_name += ".root";
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
        run.Run(fit_energies);
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
            run.Run(fit_energies);
            file->Close();
        }
    }
    return 0;
}
