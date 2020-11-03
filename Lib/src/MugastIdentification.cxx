#include "MugastIdentification.h"

#ifdef VERBOSE_DEBUG
#  define DEBUG(x, y) std::cout << (x) << (y) << std::endl;
#else
#  define DEBUG(x, y) 
#endif

MugastIdentification::MugastIdentification() :
                           light_particles({"m1_z1", "m2_z1", "m4_z2"}),
                           strips({"X", "Y"}),
                           cuts_MG({1, 3, 4, 5, 7, 11}),
                           cuts_M({1, 2, 4}),
                           cuts_Z({1, 2}),
                           fragments({"m1_z1", "m2_z1", "m4_z2", "m46_z18", "m47_z19"}),
                           layers({"ice_front",
                                   "havar_front",
                                   "he3_front",
                                   "he3_back",
                                   "havar_back",
                                   "ice_back",
                                   "al_front"}),
                           use_constant_thickness(false),
                           with_cuts(true),
                           havar_thickness(3.8E-3*UNITS::mm), //in mm
                           data(nullptr),
                           fragment(nullptr){
    system("rm -f Configs/Interpolations/TW_Ice_Thickness.root");//TODO: remove this.. it is temporary
}

MugastIdentification::~MugastIdentification()
{
    delete gas_thickness;
    delete havar_angle;
    delete data;
    delete fragment;

    for (const auto &MG : cuts_MG)
    {
        delete calibrations_TY[MG];
    }
    for (const auto &particle : light_particles)
    {
        if (particle == "m4_z2")
            continue; //TODO: generate table
        for (const auto &layer : layers)
        {
            delete energy_loss[particle][layer];
        }
    }
    for (const auto &layer : layers)
    {
        delete energy_loss["beam"][layer];
    }
}

bool MugastIdentification::Initialize(const double &beam_energy,
                                      const TVector3 &target_pos){
    DEBUG("------------>MugastIdentification::Initialize()", "");

    this->beam_energy = beam_energy;
    this->target_pos = target_pos;

    if (use_constant_thickness){
        current_ice_thickness.first = 3E-2*UNITS::mm;
        current_ice_thickness.second = 3E-2*UNITS::mm;
    }

    //Setup of reactions////////////////////////////////
    reaction.emplace("M47_Z19_m1_z1", new ReactionReconstruction2body(
                                            ReactionReconstruction2body::ReactionInput2body({
                                            ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                            ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                            ReactionFragment::FragmentSettings(1, 1, 0, 0, 0),
                                            ReactionFragment::FragmentSettings(48, 19, 0, 0, 0)})));
    reaction.emplace("M47_Z19_m2_z1", new ReactionReconstruction2body(
            ReactionReconstruction2body::ReactionInput2body({
                                                                    ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                    ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                    ReactionFragment::FragmentSettings(2, 1, 0, 0, 0),
                                                                    ReactionFragment::FragmentSettings(47, 19, 0,0, 0)})));
    reaction.emplace("M47_Z19_m4_z2", new ReactionReconstruction2body(
            ReactionReconstruction2body::ReactionInput2body({
                                                                    ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                    ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                    ReactionFragment::FragmentSettings(4, 2, 0, 0, 0),
                                                                    ReactionFragment::FragmentSettings(45, 18, 0,0, 0)})));
    beam_ref = &reaction.at("M47_Z19_m2_z1")->GetReactionFragment(1);

    //Setting up interpolations for target deformation//////
    InitializeInterpolations();

    //Cuts Initialization///////////////////////////////////
    InitializeCuts();

    //Calibration initialization////////////////////////////
    InitializeCalibration();

    //Energy Loss///////////////////////////////////////////
    InitializeELoss();

    return true;
}

bool MugastIdentification::InitializeCuts(){
    auto* tmp = new std::unordered_map<std::string, TCutG *>();
    std::vector<std::string> tmp_cut_names;

    tmp_cut_names.emplace_back("E_TOF_m1_z1_MG1");
    tmp_cut_names.emplace_back("E_TOF_m1_z1_MG3");
    tmp_cut_names.emplace_back("E_TOF_m1_z1_MG4");
    tmp_cut_names.emplace_back("E_TOF_m1_z1_MG5");
    tmp_cut_names.emplace_back("E_TOF_m1_z1_MG7");
    tmp_cut_names.emplace_back("E_TOF_m1_z1_MG11");

    tmp_cut_names.emplace_back("E_TOF_m2_z1_MG1");
    tmp_cut_names.emplace_back("E_TOF_m2_z1_MG3");
    tmp_cut_names.emplace_back("E_TOF_m2_z1_MG4");
    tmp_cut_names.emplace_back("E_TOF_m2_z1_MG5");
    tmp_cut_names.emplace_back("E_TOF_m2_z1_MG7");
    tmp_cut_names.emplace_back("E_TOF_m2_z1_MG11");

    tmp_cut_names.emplace_back("E_TOF_m4_z2_MG1");
    tmp_cut_names.emplace_back("E_TOF_m4_z2_MG3");
    tmp_cut_names.emplace_back("E_TOF_m4_z2_MG4");
    tmp_cut_names.emplace_back("E_TOF_m4_z2_MG5");
    tmp_cut_names.emplace_back("E_TOF_m4_z2_MG7");
    tmp_cut_names.emplace_back("E_TOF_m4_z2_MG11");

    for (const auto &cut_name : tmp_cut_names){
        try{
            (*tmp)[cut_name] = cuts.at(cut_name);
        }catch (const std::out_of_range &err){
            std::cerr << "Unable to find one of the E_TOF cuts : " << cut_name << "\n";
        }
    }
    cut_type["E_TOF"] = *tmp;
    return true;
}

bool MugastIdentification::InitializeCalibration(){
    for (const auto &MG : cuts_MG){
        std::string file_name = "./Configs/Calibrations/TY_MG" + to_string(MG) + ".cal";
        ifstream file(file_name);
        if (!file){
            calibrations_TY[MG] = nullptr;
            continue;
        }
        std::cout << "Calibration :" + file_name + " found\n";
        calibrations_TY[MG] = new Calibration(file_name, 1, n_strips);
    }
    return true;
}

bool MugastIdentification::InitializeELoss(){

    std::unordered_map<std::string, std::string> layer_names;
    layer_names["ice_front"]    = "Ice";
    layer_names["havar_front"]  = "Havar";
    layer_names["he3_front"]    = "Helium";
    layer_names["he3_back"]     = "Helium";
    layer_names["havar_back"]   = "Havar";
    layer_names["ice_back"]     = "Ice";
    layer_names["al_front"]     = "Al";

    std::unordered_map<std::string, std::string> fragment_names;
    fragment_names["m1_z1"]     = "Proton";
    fragment_names["m2_z1"]     = "Deuteron";
    fragment_names["m4_z2"]     = "Alpha";
    fragment_names["m46_z18"]   = "Ar46";
    fragment_names["m47_z19"]   = "K47";

    for (const auto &fragment : fragments){
        for (const auto &layer : layers){
            energy_loss[fragment][layer] = nullptr;
        }
    }

    std::string e_loss_file_path("./Configs/ELossTables/");
    std::string pressure("500Torr");

    int beam_precision = 100;
    int fragment_precision = 200;

    if (e_loss_file_path != nullptr){
        for (const auto &fragment : fragments){
            int precision;
            if (fragment == "m46_z18")
                precision = beam_precision;
            else
                precision = fragment_precision;

            for (const auto &layer : layers){
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
                                 layer_names[layer]
                                 +
                                 pressure +
                                 ".txt";
                }
                energy_loss[fragment][layer] = new EnergyLoss(tmp_string, "SRIM", precision);
            }
        }
    }

    for (const auto &layer : layers){
        energy_loss["beam"][layer] = energy_loss["m46_z18"][layer];
    }

    current_ice_thickness.first = 10E-3*UNITS::mm; //in mm
    current_ice_thickness.second = current_ice_thickness.first * ice_percentage_second;
    beam_energy_match_threashold = 0.3*UNITS::MeV; //In MeV

    std::string path_TW_vs_ice = "./Configs/Interpolations/TW_Ice_Thickness.root";
    std::ifstream tmp_file(path_TW_vs_ice);
    if (tmp_file){
        tmp_file.close();
        auto* tmp_root_file = new TFile(path_TW_vs_ice.c_str());
        TW_vs_ice_thickness = new Interpolation(tmp_root_file);
        tmp_root_file->Close();
    }else{
        TW_vs_ice_thickness = nullptr;
        TW_vs_ice.reserve(100);
    }
    return true;
}

bool MugastIdentification::InitializeInterpolations() {
    gas_thickness = new Interpolation("./Configs/Interpolations/He_thickness_500Torr.txt");
    havar_angle = new Interpolation("./Configs/Interpolations/Entrance_angle_500Torr.txt");

    auto* tmp_file = new TFile("./Configs/Interpolations/TW_Brho_M46_Z18.root");
    if (tmp_file->IsOpen()){
        TW_Brho_M46_Z18 = new Interpolation(tmp_file);
        tmp_file->Close();
        delete tmp_file;
    }else{
        std::cerr << "Unable to find TW_Brho file\n";
        TW_Brho_M46_Z18 = nullptr;
    }

    std::string tmp_file_path = "./Configs/Interpolations/TW_Ice_Thickness.root";
    std::ifstream test_file(tmp_file_path);
    if (test_file.is_open()){
        test_file.close();
        auto* tmp_file = new TFile(tmp_file_path.c_str());
        ice_thickness = new Interpolation(tmp_file);
        tmp_file->Close();
        delete tmp_file;
    }else{
        std::cerr << "Unable to find TW_Ice_Thickness file\n";
        ice_thickness = nullptr;
    }
    return true;
}

std::vector<std::pair<double, double>> MugastIdentification::GetTWvsIce(){
    return TW_vs_ice;
}

bool MugastIdentification::Identify(){
    DEBUG("------------>MugastIdentification::Identify()", "");
    //Initialization of basic structure
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii){
        fragment->Indentified[ii] = false;
        fragment->Pos[ii]   = TVector3((**(data->Mugast)).PosX[ii],
                                     (**(data->Mugast)).PosY[ii],
                                     (**(data->Mugast)).PosZ[ii]);
        fragment->TelescopeNormal[ii] = TVector3((**(data->Mugast)).TelescopeNormalX[ii],
                                                 (**(data->Mugast)).TelescopeNormalY[ii],
                                                 (**(data->Mugast)).TelescopeNormalZ[ii]);
        fragment->SI_E[ii]  = (**(data->Mugast)).DSSD_E[ii];
        fragment->SI_E2[ii] = (**(data->Mugast)).SecondLayer_E[ii];
        fragment->SI_X[ii]  = (**(data->Mugast)).DSSD_X[ii];
        fragment->SI_Y[ii]  = (**(data->Mugast)).DSSD_Y[ii];
        fragment->SI_T[ii]  = (**(data->Mugast)).DSSD_T[ii];
        fragment->T2[ii]    = (**(data->Mugast)).SecondLayer_T[ii];
        fragment->MG[ii]    = (**(data->Mugast)).TelescopeNumber[ii];
    }
    DEBUG("------------>finished: setting up fragment", "");

    //Evaluate Ice thickness
    //TODO: fix this not to make unnecessary calculations
    brho = TW_Brho_M46_Z18->Evaluate(**(data->TW));
    final_beam_energy = beam_ref->GetEkFromBrho_Q(brho, charge_state_interpolation);
    initial_beam_energy = InitialBeamEnergy(final_beam_energy, current_ice_thickness.first);

    if (TW_vs_ice_thickness == nullptr){
        if (!use_constant_thickness){
            IdentifyIceThickness();
        }
    }else{
        current_ice_thickness.first = ice_thickness->Evaluate(**(data->TW));
        current_ice_thickness.second = current_ice_thickness.first * ice_percentage_second;
    }
    for (const auto &it : reaction){
        it.second->SetBeamEnergy(MiddleTargetBeamEnergy(final_beam_energy));
    }

    DEBUG("------------>starting: recalibrations", "");
    //Applying (time) re-calibrations
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii){
        if (calibrations_TY[fragment->MG[ii]] == nullptr)
            fragment->T[ii] = fragment->SI_T[ii];
        else
            fragment->T[ii] = calibrations_TY[fragment->MG[ii]]
                                  ->Evaluate(fragment->SI_T[ii], fragment->SI_Y[ii]);
    };

    DEBUG("------------>finished: calibrations\n", "");

    //Identification with E TOF
    TCutG *tmp_cut;
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii){
        fragment->Particle[ii] = "NONE";
        for (const auto &cut_it : light_particles){
            if (with_cuts){
                try{
                    tmp_cut = cut_type["E_TOF"].at("E_TOF_" + cut_it + "_MG" +
                                                   std::to_string(static_cast<int>(fragment->MG[ii])));
                }catch (const std::out_of_range &err){
                    std::cerr << "Mugast cuts not found : "
                              << "E_TOF_" + cut_it + "_MG"
                              << std::to_string(static_cast<int>(fragment->MG[ii]))
                              << std::endl;
                    with_cuts = false;
                    continue;
                }
                if (tmp_cut->IsInside(fragment->SI_E[ii], fragment->T[ii])){
                    if (fragment->Indentified[ii])
                        throw std::runtime_error("Overlapping MUGAST E TOF gates :" +
                                                 cut_it +
                                                 "\tMG" +
                                                 fragment->MG[ii] +
                                                 "\n");

                    fragment->M[ii] = std::stoi(cut_it.substr(cut_it.find_first_of('m') + 1,
                                                              cut_it.find_first_of('_') - cut_it.find_first_of('m') - 1));

                    fragment->Z[ii] = std::stoi(cut_it.substr(cut_it.find_first_of('z') + 1,
                                                              cut_it.find_first_of('_') - cut_it.find_first_of('z') - 1));

                    fragment->Indentified[ii] = true;
                    fragment->Particle[ii] = cut_it;
                }
            }
            if (!fragment->Indentified[ii]){
                fragment->M[ii] = 0;
                fragment->Z[ii] = 0;
            }
        };
    }
    DEBUG("------------>finished: searching cuts", "");

    //Energy reconstruction
    std::unordered_map<std::string, EnergyLoss *> *ptr_tmp;
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii){
        //if (!fragment->Indentified[ii]) {
        fragment->EmissionDirection[ii] = fragment->Pos[ii] - target_pos;
        if (!fragment->Indentified[ii] || fragment->M[ii] == 4){ //TODO: fix to include alphas
            fragment->E[ii] = fragment->SI_E[ii];
            fragment->Ex[ii] = 0;
            continue;
        }
        double tmp_en = 0;

        ptr_tmp = &energy_loss["m" +
                               std::to_string(fragment->M[ii]) +
                               "_z" +
                               std::to_string(fragment->Z[ii])];

        double theta = fragment->EmissionDirection[ii].Angle(TVector3(0, 0, -1));

        //Passivation layer
        tmp_en = (*ptr_tmp)["al_front"]
                     ->EvaluateInitialEnergy(fragment->SI_E[ii],
                                             0.4E-3*UNITS::mm, //Units in mm!
                                             fragment->EmissionDirection[ii]
                                                 .Angle(fragment->TelescopeNormal[ii]));
        tmp_en = (*ptr_tmp)["ice_front"]
                     ->EvaluateInitialEnergy(tmp_en,
                                             current_ice_thickness.first,
                                             havar_angle->Evaluate(theta));
        tmp_en = (*ptr_tmp)["havar_front"]
                     ->EvaluateInitialEnergy(tmp_en,
                                             havar_thickness,
                                             havar_angle->Evaluate(theta));
        tmp_en = (*ptr_tmp)["he3_front"]
                     ->EvaluateInitialEnergy(tmp_en,
                                             gas_thickness->Evaluate(theta)*UNITS::mm,
                                             0.);

        fragment->E[ii] = tmp_en;

        if ((reaction_it = reaction.find("M" + std::to_string(data->VAMOS_id_M) +
                                         "_Z" + std::to_string(data->VAMOS_id_Z) +
                                         "_" +
                                         fragment->Particle[ii])) != reaction.end()){
            fragment->Ex[ii] = reaction_it->second->Set_Ek_Theta(fragment->E[ii], fragment->EmissionDirection[ii].Theta());//WARNING this could be wring
            fragment->E_CM[ii] = reaction_it->second->GetReactionFragment(3).Get_Ek_cm();
        }
        else{
            fragment->Ex[ii] = 0;
        }
    }

    DEBUG("------------>Finished : MugastIdentification::Identify()", "");

    return true;
}

void MugastIdentification::IdentifyIceThickness(){
    DEBUG("------------>Started: MugastIdentification::IdentifyIceThickness()", "");
    DEBUG("Final beam energy : " , final_beam_energy);
    DEBUG("Computed initial beam energy : " , initial_beam_energy);
    DEBUG(" With thickness : " , current_ice_thickness.first);
    DEBUG(" With brho : " , brho);
    DEBUG(" and mass : " , beam_ref->Get_M2());

    if (abs(beam_energy - initial_beam_energy) > beam_energy_match_threashold){
        DEBUG("Before minimizer call, ice_thickness.first = ", current_ice_thickness.first);
        double tmp_threashold{2E-4*UNITS::mm};
        double tmp_precision{0.1*UNITS::MeV};

        ice_thickness_minimizer = new Minimizer([this](const double &tck) { return pow(this->beam_energy - this->InitialBeamEnergy(this->final_beam_energy, tck), 2)/(this->beam_energy*this->beam_energy); },
                                                current_ice_thickness.first,        //starting value
                                                3E-4 * current_ice_thickness.first, //learning rate
                                                tmp_threashold,                     //threashold
                                                100,                                //max_steps
                                                1,                                  //quenching
                                                1.E-4*current_ice_thickness.first);//h

        while (fabs(beam_energy - InitialBeamEnergy(final_beam_energy, current_ice_thickness.first)) > tmp_precision){
            tmp_threashold/=2;
            ice_thickness_minimizer->SetThreashold(tmp_threashold);
            current_ice_thickness.first = ice_thickness_minimizer->Minimize();
            current_ice_thickness.second = current_ice_thickness.first * ice_percentage_second;
        }

        DEBUG("After minimizer call, ice_thickness.first = ", current_ice_thickness.first);
        DEBUG(" energy difference : ", beam_energy - InitialBeamEnergy(final_beam_energy, current_ice_thickness.first));
    }
    DEBUG("------------>Finished: MugastIdentification::IdentifyIceThickness()", "");
}

double MugastIdentification::InitialBeamEnergy(double beam_energy_from_brho, double tck){
    DEBUG("------------>Started: MugastIdentification::InitialBeamEnergy(), thickness : ", tck);
    beam_energy_from_brho =
        energy_loss["beam"]["ice_back"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    tck * ice_percentage_second,
                                    0);

    beam_energy_from_brho =
        energy_loss["beam"]["havar_back"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    havar_thickness,
                                    0);

    beam_energy_from_brho =
        energy_loss["beam"]["he3_front"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    average_beam_thickness,
                                    0);

    beam_energy_from_brho =
        energy_loss["beam"]["havar_front"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    havar_thickness,
                                    0);

    beam_energy_from_brho =
        energy_loss["beam"]["ice_front"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    tck,
                                    0);
    DEBUG("------------>Finished: MugastIdentification::InitialBeamEnergy()", "");
    return beam_energy_from_brho;
}

double MugastIdentification::MiddleTargetBeamEnergy(double beam_energy_from_brho){
    beam_energy_from_brho =
        energy_loss["beam"]["ice_back"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    current_ice_thickness.first,
                                    0);

    beam_energy_from_brho =
        energy_loss["beam"]["havar_back"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    havar_thickness,
                                    0);

    beam_energy_from_brho =
        energy_loss["beam"]["he3_front"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    average_beam_thickness / 2.,
                                    0);

    return beam_energy_from_brho;
}