#include "MugastIdentification.h"

MugastIdentification::MugastIdentification() : cuts_MG({1, 3, 4, 5, 7, 11}),
                                               cuts_M({1, 2, 4}),
                                               cuts_Z({1, 2}) {}

MugastIdentification::~MugastIdentification() {
    delete gas_thickness;
    delete havar_angle;
    delete data;
    delete fragment;
}

bool MugastIdentification::Initialize() {
    gas_thickness = new Interpolation("./Configs/Interpolations/GasThickness.txt");
    havar_angle = new Interpolation("./Configs/Interpolations/EntranceAngleHavar.txt");

    //Masses in MeV [M][Z] as ints
    for (const auto &it_M : cuts_M) {
        for (const auto &it_Z : cuts_Z) {
            mass[it_M][it_Z] = it_M * AMU_TO_MEV;
        }
    }

    //More precise masses [M][Z]
    mass[2][1] = 2.01410177812 * AMU_TO_MEV;  //In MeV
    mass[1][1] = 1.00782503223 * AMU_TO_MEV;  //In MeV

    std::unordered_map<std::string, TCutG *> *tmp 
        = new std::unordered_map<std::string, TCutG *>();

    try {
        //Protons
        (*tmp)["E_TOF_m1_z1_MG1"]   = cuts.at("E_TOF_m1_z1_MG1");
        (*tmp)["E_TOF_m1_z1_MG3"]   = cuts.at("E_TOF_m1_z1_MG3");
        (*tmp)["E_TOF_m1_z1_MG4"]   = cuts.at("E_TOF_m1_z1_MG4");
        (*tmp)["E_TOF_m1_z1_MG5"]   = cuts.at("E_TOF_m1_z1_MG5");
        (*tmp)["E_TOF_m1_z1_MG7"]   = cuts.at("E_TOF_m1_z1_MG7");
        (*tmp)["E_TOF_m1_z1_MG11"]  = cuts.at("E_TOF_m1_z1_MG11");
        //Deuterons
        (*tmp)["E_TOF_m2_z1_MG1"]   = cuts.at("E_TOF_m2_z1_MG1");
        (*tmp)["E_TOF_m2_z1_MG3"]   = cuts.at("E_TOF_m2_z1_MG3");
        (*tmp)["E_TOF_m2_z1_MG4"]   = cuts.at("E_TOF_m2_z1_MG4");
        (*tmp)["E_TOF_m2_z1_MG5"]   = cuts.at("E_TOF_m2_z1_MG5");
        (*tmp)["E_TOF_m2_z1_MG7"]   = cuts.at("E_TOF_m2_z1_MG7");
        (*tmp)["E_TOF_m2_z1_MG11"]  = cuts.at("E_TOF_m2_z1_MG11");
    } catch (const std::out_of_range &err) {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the E_TOF cuts\n");
    }

    cut_type["E_TOF"] = *tmp;

    for (const auto &cut_type_it : cut_type) {
        for (const auto &cut_it : cut_type_it) {
            cut_detector[std::stoi(cut_it.first.substr(cut_it.first.find_last_of("_") + 3))] = cut_it.second;
        }
    }

    return true;
}