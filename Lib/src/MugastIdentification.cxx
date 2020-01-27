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

bool MugastIdentification::Initialize(const double & beam_energy,
                                        const TVector3 & target_pos) {
    gas_thickness = new Interpolation("./Configs/Interpolations/GasThickness.txt");
    havar_angle = new Interpolation("./Configs/Interpolations/EntranceAngleHavar.txt");

    this->beam_energy = beam_energy;
    this->target_pos = target_pos;
    //Cuts Initialization///////////////////////////////////
    InitializeCuts();

    //Calibration initialization//////////////////////////////////////////
    InitializeCalibration();

    //Energy Loss//////////////////////////////////////////////////////////
    InitializeELoss();

    return true;
}

bool MugastIdentification::InitializeCuts() {
    //Masses in MeV [M][Z] as ints
    for (const auto &it_M : cuts_M) {
        for (const auto &it_Z : cuts_Z) {
            mass[it_M][it_Z] = it_M * AMU_TO_MEV;
        }
    }

    //More precise masses [M][Z]
    mass[2][1] = 2.01410177812 * AMU_TO_MEV;  //In MeV
    mass[1][1] = 1.00782503223 * AMU_TO_MEV;  //In MeV

    std::unordered_map<std::string, TCutG *> *tmp = new std::unordered_map<std::string, TCutG *>();

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

    for (const auto &cut_name : tmp_cut_names) {
        try {
            (*tmp)[cut_name] = cuts.at(cut_name);
        } catch (const std::out_of_range &err) {
            std::cerr << "Unable to find one of the E_TOF cuts : " << cut_name << "\n";
        }
    }

    cut_type["E_TOF"] = *tmp;

    return true;
}

bool MugastIdentification::InitializeCalibration() {
    for (const auto &MG : cuts_MG) {
        std::string file_name = "./Configs/Calibrations/TY_MG" + to_string(MG) + ".cal";
        ifstream file(file_name);
        if (!file) {
            calibrations_TY[MG] = nullptr;
            continue;
        }
        std::cout << "Calibration :" + file_name + " found\n";
        calibrations_TY[MG] = new Calibration(file_name, 1, n_strips);
    }

    return true;
}

bool MugastIdentification::InitializeELoss() {
    std::vector<std::string> layers = {"ice_front",
                                       "havar_front",
                                       "he3_front",
                                       "he3_back",
                                       "havar_back",
                                       "ice_back",
                                       "al_front"};

    std::unordered_map<std::string, std::string> layer_names;

    layer_names["ice_front"] = "Ice";
    layer_names["havar_front"] = "Havar";
    layer_names["he3_front"] = "He_gas";
    layer_names["he3_back"] = "Ice";
    layer_names["havar_back"] = "Havar";
    layer_names["ice_back"] = "Ice";
    layer_names["al_front"] = "Al";

    std::unordered_map<std::string, std::string> particle_names;

    particle_names["m1_z1"] = "proton";
    particle_names["m2_z1"] = "deuteron";
    particle_names["m4_z2"] = "alpha";

    for (const auto &particle : particles) {
        for (const auto &layer : layers) {
            energy_loss[particle][layer] = nullptr;
        }
    }

    std::ifstream configs_file("./Configs/Configs.txt");
    std::string line;

    std::string tmp_string;

    std::string e_loss_file_path("not_found");

    if (configs_file) {
        while (std::getline(configs_file, line)) {
            std::istringstream line_str(line);
            line_str >> tmp_string;
            if (tmp_string.compare("ELossTableFolder:") == 0) {
                line_str >> e_loss_file_path;
                break;
            }
        }
    }

    if (e_loss_file_path != 0) {
        for (const auto &particle : particles) {
            if (particle == "m4_z2") continue;  //TODO: generate table
            for (const auto &layer : layers) {
                tmp_string = e_loss_file_path +
                             particle_names[particle] +
                             "_" +
                             layer_names[layer] +
                             ".G4table";
                energy_loss[particle][layer] =
                    new NPL::EnergyLoss(tmp_string, "G4Table", 1000);
            }
        }
    }

    for (const auto &layer : layers) {
        energy_loss["beam"][layer] =
            new NPL::EnergyLoss("46Ar_" + layer_names[layer] + "G4table",
                                "G4Table",
                                1000);
    }
    return true;
}