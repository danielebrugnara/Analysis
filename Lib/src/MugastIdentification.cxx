#include "MugastIdentification.h"

MugastIdentification::MugastIdentification() : cuts_MG({1, 3, 4, 5, 7, 11}),
                                               cuts_M({1, 2, 4}),
                                               cuts_Z({1, 2}),
                                               particles({"m1_z1", "m2_z1", "m4_z2"}),
                                               strips({"X", "Y"}),
                                               data(nullptr),
                                               fragment(nullptr),
                                               with_cuts(true) {}

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

    std::vector<std::string> tmp_cut_names;
    tmp_cut_names.push_back("E_TOF_m1_z1_MG1");
    tmp_cut_names.push_back("E_TOF_m1_z1_MG3");
    tmp_cut_names.push_back("E_TOF_m1_z1_MG4");
    tmp_cut_names.push_back("E_TOF_m1_z1_MG5");
    tmp_cut_names.push_back("E_TOF_m1_z1_MG7");
    tmp_cut_names.push_back("E_TOF_m1_z1_MG11");
                            
    tmp_cut_names.push_back("E_TOF_m2_z1_MG1");
    tmp_cut_names.push_back("E_TOF_m2_z1_MG3");
    tmp_cut_names.push_back("E_TOF_m2_z1_MG4");
    tmp_cut_names.push_back("E_TOF_m2_z1_MG5");
    tmp_cut_names.push_back("E_TOF_m2_z1_MG7");
    tmp_cut_names.push_back("E_TOF_m2_z1_MG11");

    for (const auto & cut_name : tmp_cut_names){
        try{
            (*tmp)[cut_name] = cuts.at(cut_name);
        }catch (const std::out_of_range &err) {
            std::cerr << "Unable to find one of the E_TOF cuts : "<< cut_name << "\n";
        }
    }

    cut_type["E_TOF"] = *tmp;

    return true;
}
