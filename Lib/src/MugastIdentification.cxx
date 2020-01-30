#include "MugastIdentification.h"

MugastIdentification::MugastIdentification() : cuts_MG({1, 3, 4, 5, 7, 11}),
                                               cuts_M({1, 2, 4}),
                                               cuts_Z({1, 2}),
                                               particles({"m1_z1", "m2_z1", "m4_z2"}),
                                               strips({"X", "Y"}),
                                               layers({"ice_front",
                                                       "havar_front",
                                                       "he3_front",
                                                       "he3_back",
                                                       "havar_back",
                                                       "ice_back",
                                                       "al_front"}), 
                                               data(nullptr),
                                               fragment(nullptr),
                                               with_cuts(true),
                                               havar_thickness(3.8E-3){//in mm
}

MugastIdentification::~MugastIdentification() {
    delete gas_thickness;
    delete havar_angle;
    delete data;
    delete fragment;
    //delete reaction;
    for (const auto &MG : cuts_MG) {
        delete calibrations_TY[MG];
    }
    for (const auto &particle : particles) {
        if (particle == "m4_z2") continue;  //TODO: generate table
        for (const auto &layer : layers) {
            delete energy_loss[particle][layer];
        }
    }
    for (const auto &layer : layers) {
        delete energy_loss["beam"][layer];
    }
}

bool MugastIdentification::Initialize(const double &beam_energy,
                                      const TVector3 &target_pos) {
#ifdef VERBOSE_DEBUG
    std::cout << "------------>MugastIdentification::Initialize()\n";
#endif
    reaction["M47_Z19_m1_z1"] = new NPL::Reaction("46Ar(3He,p)48K@" + std::to_string(beam_energy));
    reaction["M47_Z19_m2_z1"] = new NPL::Reaction("46Ar(3He,d)47K@" + std::to_string(beam_energy));
    reaction["M47_Z19_m4_z2"] = new NPL::Reaction("46Ar(3He,4He)45Ar@" + std::to_string(beam_energy));


    gas_thickness = new Interpolation("./Configs/Interpolations/GasThickness.txt");
    havar_angle = new Interpolation("./Configs/Interpolations/EntranceAngleHavar.txt");

    TFile *tmp_file = new TFile("./Configs/Interpolations/TW_Brho_M46_Z18.root");
    if (tmp_file){
        TW_Brho_M46_Z18 = new Interpolation(tmp_file);
        tmp_file->Close();
        delete tmp_file;
    }else{
        TW_Brho_M46_Z18 = nullptr;
    }

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
    mass[46][18] = 45.968082712 * AMU_TO_MEV;  //in MeV

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

    tmp_cut_names.push_back("E_TOF_m4_z2_MG1");
    tmp_cut_names.push_back("E_TOF_m4_z2_MG3");
    tmp_cut_names.push_back("E_TOF_m4_z2_MG4");
    tmp_cut_names.push_back("E_TOF_m4_z2_MG5");
    tmp_cut_names.push_back("E_TOF_m4_z2_MG7");
    tmp_cut_names.push_back("E_TOF_m4_z2_MG11");

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
            new NPL::EnergyLoss(e_loss_file_path+"Ar46_" + 
                                    layer_names[layer] + 
                                    ".G4table",
                                "G4Table",
                                1000);
    }

    energy_loss["beam"]["ice_front"] =
        //new NPL::EnergyLoss("/home/daniele/Projects/Analysis/Configs/ELossTables/Argon_in_ice.txt", "SRIM", 1000);
        new NPL::EnergyLoss("./Configs/ELossTables/Argon_in_ice.txt", "SRIM", 1000);
    energy_loss["beam"]["ice_back"] =
        //new NPL::EnergyLoss("/home/daniele/Projects/Analysis/Configs/ELossTables/Argon_in_ice.txt", "SRIM", 1000);
        new NPL::EnergyLoss("./Configs/ELossTables/Argon_in_ice.txt", "SRIM", 1000);

    current_ice_thickness = 10E-3;
    beam_energy_match_threashold = 0.3; //In MeV
    return true;
}
