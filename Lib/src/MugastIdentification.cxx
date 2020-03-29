#include "MugastIdentification.h"

MugastIdentification::MugastIdentification() : cuts_MG({1, 3, 4, 5, 7, 11}),
                                               cuts_M({1, 2, 4}),
                                               cuts_Z({1, 2}),
                                               light_particles({"m1_z1", "m2_z1", "m4_z2"}),
                                               fragments({"m1_z1", "m2_z1", "m4_z2","m46_z18", "m47_z19"}),
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
    for (const auto &particle : light_particles) {
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


    gas_thickness = new Interpolation("./Configs/Interpolations/He_thickness_500Torr.txt");
    havar_angle = new Interpolation("./Configs/Interpolations/Entrance_angle_500Torr.txt");

    TFile *tmp_file = new TFile("./Configs/Interpolations/TW_Brho_M46_Z18.root");
    if (tmp_file){
        TW_Brho_M46_Z18 = new Interpolation(tmp_file);
        tmp_file->Close();
        delete tmp_file;
    }else{
        TW_Brho_M46_Z18 = nullptr;
    }

    std::string tmp_file_path = "./Configs/Interpolations/TW_Ice_Thickness.root";
    std::ifstream test_file(tmp_file_path);
    if (test_file){
        test_file.close();
        TFile* tmp_file = new TFile(tmp_file_path.c_str());
        ice_thickness = new Interpolation(tmp_file);
        tmp_file->Close();
        delete tmp_file;
    }else{
        ice_thickness = nullptr;
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
    layer_names["he3_front"] = "Helium";
    layer_names["he3_back"] = "Helium";
    layer_names["havar_back"] = "Havar";
    layer_names["ice_back"] = "Ice";
    layer_names["al_front"] = "Al";

    std::unordered_map<std::string, std::string> fragment_names;

    fragment_names["m1_z1"] = "Proton";
    fragment_names["m2_z1"] = "Deuteron";
    fragment_names["m4_z2"] = "Alpha";
    fragment_names["m46_z18"] = "Ar46";
    fragment_names["m47_z19"] = "K47";

    for (const auto &fragment : fragments) {
        for (const auto &layer : layers) {
            energy_loss[fragment][layer] = nullptr;
        }
    }

    std::string e_loss_file_path("./Configs/ELossTables/");
    std::string pressure("500Torr");

    int beam_precision = 100;
    int fragment_precision = 200;

    
    if (e_loss_file_path != 0) {
        for (const auto &fragment : fragments) {

            int precision;
            if (fragment == "m46_z18")
                precision = beam_precision;
            else
                precision = fragment_precision;
            
            for (const auto &layer : layers) {
                std::string tmp_string;
                if (layer_names[layer] != "Helium"){
                    tmp_string = e_loss_file_path +
                                                fragment_names[fragment] +
                                                "_in_" +
                                                layer_names[layer] +
                                                ".txt";
                }else{
                    tmp_string = e_loss_file_path +
                                                fragment_names[fragment] +
                                                "_in_" +
                                                layer_names[layer] +
                                                pressure+
                                                ".txt";
                }
                energy_loss[fragment][layer] =
                    new EnergyLoss(tmp_string, "SRIM", precision);
            }
        }
    }

    for (const auto &layer : layers) {
        energy_loss["beam"][layer] = energy_loss["m46_z18"][layer];
    }

    current_ice_thickness.first = 10E-3;
    current_ice_thickness.second = current_ice_thickness.first * ice_percentage_second;
    beam_energy_match_threashold = 0.3; //In MeV
    
    std::string path_TW_vs_ice = "./Configs/Interpolations/TW_Ice_Thickness.root";
    std::ifstream tmp_file(path_TW_vs_ice);
    if (tmp_file){
        tmp_file.close();
        TFile *tmp_root_file = new TFile(path_TW_vs_ice.c_str());
        TW_vs_ice_thickness = new Interpolation(tmp_root_file);
        tmp_root_file->Close();
    }else{
        TW_vs_ice_thickness = nullptr;
        TW_vs_ice.reserve(100);
    }
    return true;
}

std::vector<std::pair<double, double>> MugastIdentification::GetTWvsIce(){
    return TW_vs_ice;
}
