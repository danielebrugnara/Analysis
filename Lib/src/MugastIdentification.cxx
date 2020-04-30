#include "MugastIdentification.h"

MugastIdentification::MugastIdentification() : cuts_MG({1, 3, 4, 5, 7, 11}),
                                               cuts_M({1, 2, 4}),
                                               cuts_Z({1, 2}),
                                               light_particles({"m1_z1", "m2_z1", "m4_z2"}),
                                               fragments({"m1_z1", "m2_z1", "m4_z2", "m46_z18", "m47_z19"}),
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
                                               use_constant_thickness(false),
                                               with_cuts(true),
                                               havar_thickness(3.8E-3)
{ //in mm
}

MugastIdentification::~MugastIdentification()
{
    delete gas_thickness;
    delete havar_angle;
    delete data;
    delete fragment;
    //delete reaction;
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
                                      const TVector3 &target_pos)
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>MugastIdentification::Initialize()\n";
#endif

    if (use_constant_thickness)
    {
        current_ice_thickness.first = 3E-2;
        current_ice_thickness.second = 3E-2;
    }

    reaction["M47_Z19_m1_z1"] = new NPL::Reaction("46Ar(3He,p)48K@" + std::to_string(beam_energy));
    reaction["M47_Z19_m2_z1"] = new NPL::Reaction("46Ar(3He,d)47K@" + std::to_string(beam_energy));
    reaction["M47_Z19_m4_z2"] = new NPL::Reaction("46Ar(3He,4He)45Ar@" + std::to_string(beam_energy));

    //Masses in MeV [M][Z] as ints
    for (const auto &it_M : cuts_M)
    {
        for (const auto &it_Z : cuts_Z)
        {
            mass[it_M][it_Z] = NPL::Nucleus(it_Z, it_M).Mass();
            std::cout << "MASS::::::::: " << mass[it_M][it_Z] << std::endl;
        }
    }
    mass[46][18] = NPL::Nucleus(18, 46).Mass();

    gas_thickness = new Interpolation("./Configs/Interpolations/He_thickness_500Torr.txt");
    havar_angle = new Interpolation("./Configs/Interpolations/Entrance_angle_500Torr.txt");

    TFile *tmp_file = new TFile("./Configs/Interpolations/TW_Brho_M46_Z18.root");
    if (tmp_file)
    {
        TW_Brho_M46_Z18 = new Interpolation(tmp_file);
        tmp_file->Close();
        delete tmp_file;
    }
    else
    {
        TW_Brho_M46_Z18 = nullptr;
    }

    std::string tmp_file_path = "./Configs/Interpolations/TW_Ice_Thickness.root";
    std::ifstream test_file(tmp_file_path);
    if (test_file)
    {
        test_file.close();
        TFile *tmp_file = new TFile(tmp_file_path.c_str());
        ice_thickness = new Interpolation(tmp_file);
        tmp_file->Close();
        delete tmp_file;
    }
    else
    {
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

bool MugastIdentification::InitializeCuts()
{
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

    for (const auto &cut_name : tmp_cut_names)
    {
        try
        {
            (*tmp)[cut_name] = cuts.at(cut_name);
        }
        catch (const std::out_of_range &err)
        {
            std::cerr << "Unable to find one of the E_TOF cuts : " << cut_name << "\n";
        }
    }

    cut_type["E_TOF"] = *tmp;

    return true;
}

bool MugastIdentification::InitializeCalibration()
{
    for (const auto &MG : cuts_MG)
    {
        std::string file_name = "./Configs/Calibrations/TY_MG" + to_string(MG) + ".cal";
        ifstream file(file_name);
        if (!file)
        {
            calibrations_TY[MG] = nullptr;
            continue;
        }
        std::cout << "Calibration :" + file_name + " found\n";
        calibrations_TY[MG] = new Calibration(file_name, 1, n_strips);
    }

    return true;
}

bool MugastIdentification::InitializeELoss()
{
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

    for (const auto &fragment : fragments)
    {
        for (const auto &layer : layers)
        {
            energy_loss[fragment][layer] = nullptr;
        }
    }

    std::string e_loss_file_path("./Configs/ELossTables/");
    std::string pressure("500Torr");

    int beam_precision = 100;
    int fragment_precision = 200;

    if (e_loss_file_path != 0)
    {
        for (const auto &fragment : fragments)
        {

            int precision;
            if (fragment == "m46_z18")
                precision = beam_precision;
            else
                precision = fragment_precision;

            for (const auto &layer : layers)
            {
                std::string tmp_string;
                if (layer_names[layer] != "Helium")
                {
                    tmp_string = e_loss_file_path +
                                 fragment_names[fragment] +
                                 "_in_" +
                                 layer_names[layer] +
                                 ".txt";
                }
                else
                {
                    tmp_string = e_loss_file_path +
                                 fragment_names[fragment] +
                                 "_in_" +
                                 layer_names[layer] +
                                 pressure +
                                 ".txt";
                }
                energy_loss[fragment][layer] =
                    new EnergyLoss(tmp_string, "SRIM", precision);
            }
        }
    }

    for (const auto &layer : layers)
    {
        energy_loss["beam"][layer] = energy_loss["m46_z18"][layer];
    }

    current_ice_thickness.first = 10E-3; //in mm
    current_ice_thickness.second = current_ice_thickness.first * ice_percentage_second;
    beam_energy_match_threashold = 0.3; //In MeV

    std::string path_TW_vs_ice = "./Configs/Interpolations/TW_Ice_Thickness.root";
    std::ifstream tmp_file(path_TW_vs_ice);
    if (tmp_file)
    {
        tmp_file.close();
        TFile *tmp_root_file = new TFile(path_TW_vs_ice.c_str());
        TW_vs_ice_thickness = new Interpolation(tmp_root_file);
        tmp_root_file->Close();
    }
    else
    {
        TW_vs_ice_thickness = nullptr;
        TW_vs_ice.reserve(100);
    }
    return true;
}

std::vector<std::pair<double, double>> MugastIdentification::GetTWvsIce()
{
    return TW_vs_ice;
}

bool MugastIdentification::Identify()
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>MugastIdentification::Identify()\n";
#endif
    //Initialization of basic structure
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii)
    {
        fragment->Indentified[ii] = false;
        fragment->Pos[ii] = TVector3((**(data->Mugast)).PosX[ii],
                                     (**(data->Mugast)).PosY[ii],
                                     (**(data->Mugast)).PosZ[ii]);
        fragment->TelescopeNormal[ii] = TVector3((**(data->Mugast)).TelescopeNormalX[ii],
                                                 (**(data->Mugast)).TelescopeNormalY[ii],
                                                 (**(data->Mugast)).TelescopeNormalZ[ii]);
        fragment->SI_E[ii] = (**(data->Mugast)).DSSD_E[ii];
        fragment->SI_E2[ii] = (**(data->Mugast)).SecondLayer_E[ii];
        fragment->SI_X[ii] = (**(data->Mugast)).DSSD_X[ii];
        fragment->SI_Y[ii] = (**(data->Mugast)).DSSD_Y[ii];
        fragment->SI_T[ii] = (**(data->Mugast)).DSSD_T[ii];
        fragment->T2[ii] = (**(data->Mugast)).SecondLayer_T[ii];
        fragment->MG[ii] = (**(data->Mugast)).TelescopeNumber[ii];
    }
#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: setting up fragment\n";
#endif

    //Evaluate Ice thickness
    if (TW_vs_ice_thickness == nullptr)
    {
        if (use_constant_thickness)
        {
        }
        else
        {
            IdentifyIceThickness();
        }
    }
    else
    {
        current_ice_thickness.first = ice_thickness->Evaluate(**(data->TW));
        current_ice_thickness.second = current_ice_thickness.first * ice_percentage_second;
    }

    //Applying (time) re-calibrations
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii)
    {
        if (calibrations_TY[fragment->MG[ii]] == nullptr)
            fragment->T[ii] = fragment->SI_T[ii];
        else
            fragment->T[ii] = calibrations_TY[fragment->MG[ii]]
                                  ->Evaluate(fragment->SI_T[ii], fragment->SI_Y[ii]);
    };
#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: calibrations\n";
#endif
    //Identification with E TOF
    TCutG *tmp_cut;
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii)
    {
        fragment->Particle[ii] = "NONE";
        for (const auto &cut_it : light_particles)
        {
            if (with_cuts)
            {
                try
                {
                    tmp_cut = cut_type["E_TOF"].at("E_TOF_" + cut_it + "_MG" +
                                                   std::to_string(static_cast<int>(fragment->MG[ii])));
                }
                catch (const std::out_of_range &err)
                {
                    std::cerr << "Mugast cuts not found : "
                              << "E_TOF_" + cut_it + "_MG"
                              << std::to_string(static_cast<int>(fragment->MG[ii]))
                              << std::endl;
                    with_cuts = false;
                    continue;
                }
                if (tmp_cut->IsInside(fragment->SI_E[ii], fragment->T[ii]))
                {
                    if (fragment->Indentified[ii])
                        throw std::runtime_error("Overlapping MUGAST E TOF gates :" +
                                                 cut_it +
                                                 "\tMG" +
                                                 fragment->MG[ii] +
                                                 "\n");

                    fragment->M[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("m") + 1,
                                                              cut_it.find_first_of("_") - cut_it.find_first_of("m") - 1));

                    fragment->Z[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("z") + 1,
                                                              cut_it.find_first_of("_") - cut_it.find_first_of("z") - 1));

                    fragment->Indentified[ii] = true;
                    fragment->Particle[ii] = cut_it;
                }
            }
            if (!fragment->Indentified[ii])
            {
                fragment->M[ii] = 0;
                fragment->Z[ii] = 0;
            }
        };
    }
#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: searching cuts\n";
#endif

    //Energy reconstruction
    std::unordered_map<std::string, EnergyLoss *> *ptr_tmp;
    for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii)
    {
        //if (!fragment->Indentified[ii]) {
        fragment->EmissionDirection[ii] = fragment->Pos[ii] - target_pos;
        if (!fragment->Indentified[ii] || fragment->M[ii] == 4)
        { //TODO: fix to include alphas
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
                                             0.4E-3, //Units in mm!
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
                                             gas_thickness->Evaluate(theta),
                                             0.);

        fragment->E[ii] = tmp_en;

        if ((reaction_it = reaction.find("M" + std::to_string(data->VAMOS_id_M) +
                                         "_Z" + std::to_string(data->VAMOS_id_Z) +
                                         "_" +
                                         fragment->Particle[ii])) != reaction.end())
        {
            fragment->Ex[ii] = reaction_it->second
                                   ->ReconstructRelativistic(fragment->E[ii],
                                                             Get_ThetaLab(ii));
        }
        else
        {
            fragment->Ex[ii] = 0;
        }
    }

    //Ex

#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished : MugastIdentification::Identify()\n";
#endif

    //TODO: compute Ex here
    return true;
}

void MugastIdentification::IdentifyIceThickness()
{

#ifdef VERBOSE_DEBUG
    std::cout << "------------>Started: MugastIdentification::IdentifyIceThickness()\n";
#endif

    brho = TW_Brho_M46_Z18->Evaluate(**(data->TW));
    final_beam_energy = sqrt(pow(brho / 3.3356E-3 * charge_state_interpolation, 2) + pow(mass[46][18], 2)) - (mass[46][18]);

    initial_beam_energy = InitialBeamEnergy(final_beam_energy, current_ice_thickness.first);

#ifdef VERBOSE_DEBUG
    std::cout << "Final beam energy : " << final_beam_energy << "\n";
    std::cout << "Computed initial beam energy : " << initial_beam_energy
              << " With thickness : " << current_ice_thickness.first
              << " With brho : " << brho
              << " and mass : " << mass[46][18] << std::endl;
#endif

    if (abs(beam_energy - initial_beam_energy) > beam_energy_match_threashold)
    {
        //ice_thickness_minimizer = new Minimizer(this->Thickness_discr,
        //                                        current_ice_thickness.first, 0.7, 0.1, 100, 1, 1E-3);
#ifdef VERBOSE_DEBUG
        std::cout << "Before minimizer call, ice_thickness.first = " << current_ice_thickness.first << "\n";
#endif
        double tmp_threashold{2E-4};
        double tmp_precision{0.1}; //in MeV                                              //MeV
        ice_thickness_minimizer = new Minimizer([this](const double &tck) { return pow(this->beam_energy - this->InitialBeamEnergy(this->final_beam_energy, tck), 2); },
                                                current_ice_thickness.first,        //starting value
                                                1E-6 * current_ice_thickness.first, //learning rate
                                                tmp_threashold, //threashold
                                                100,                                //max_steps
                                                1,                                  //quenching
                                                1E-4);                              //h
        while (fabs(beam_energy - InitialBeamEnergy(final_beam_energy, current_ice_thickness.first)) > tmp_precision)
        {
            tmp_threashold/=2;
            ice_thickness_minimizer->SetThreashold(tmp_threashold);
            current_ice_thickness.first = ice_thickness_minimizer->Minimize();
            current_ice_thickness.second = current_ice_thickness.first * ice_percentage_second;
        }

        //#ifdef VERBOSE_DEBUG
        std::cout << "After minimizer call, ice_thickness.first = " << current_ice_thickness.first
                  << " energy difference : " << beam_energy - InitialBeamEnergy(final_beam_energy, current_ice_thickness.first)
                  << "\n";
        //#endif

        //        ice_thickness_minimizer = new Minimizer((beam_energy - initial_beam_energy) / gradient_descent_normalization,
        //                                                current_ice_thickness.first,
        //                                                7E-7,
        //                                                0.002,
        //                                                beam_energy_match_threashold / gradient_descent_normalization,
        //                                                200,
        //                                                1);
        //        double energy_difference;
        //        for (int ii = 0; ii < 100; ++ii)
        //        {
        //            energy_difference = beam_energy - InitialBeamEnergy(final_beam_energy);
        //            if (abs(energy_difference) < beam_energy_match_threashold)
        //                break;
        //            current_ice_thickness.first =
        //                ice_thickness_minimizer
        //                    ->PerformStep(energy_difference / gradient_descent_normalization);
        //            current_ice_thickness.second = ice_percentage_second * current_ice_thickness.first;
        //        }
        for (const auto &it : reaction)
        {
            it.second->SetBeamEnergy(MiddleTargetBeamEnergy(final_beam_energy));
        }
    }
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: MugastIdentification::IdentifyIceThickness()\n";
#endif
}

double MugastIdentification::InitialBeamEnergy(double beam_energy_from_brho, double tck)
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Started: MugastIdentification::InitialBeamEnergy(), thickness : " << tck << " \n";
#endif
    beam_energy_from_brho =
        energy_loss["beam"]["ice_back"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    tck * ice_percentage_second,
                                    0);
    //std::cout << beam_energy_from_brho << std::endl;

    beam_energy_from_brho =
        energy_loss["beam"]["havar_back"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    havar_thickness,
                                    0);
    //std::cout << beam_energy_from_brho << std::endl;

    beam_energy_from_brho =
        energy_loss["beam"]["he3_front"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    average_beam_thickness,
                                    0);
    //std::cout << beam_energy_from_brho << std::endl;

    beam_energy_from_brho =
        energy_loss["beam"]["havar_front"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    havar_thickness,
                                    0);
    //std::cout << beam_energy_from_brho << std::endl;

    beam_energy_from_brho =
        energy_loss["beam"]["ice_front"]
            ->EvaluateInitialEnergy(beam_energy_from_brho,
                                    tck,
                                    0);
    //std::cout << beam_energy_from_brho << std::endl;
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: MugastIdentification::InitialBeamEnergy()\n";
#endif
    return beam_energy_from_brho;
}

double MugastIdentification::MiddleTargetBeamEnergy(double beam_energy_from_brho)
{
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