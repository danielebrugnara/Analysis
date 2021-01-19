#include "MugastIdentification.h"

#ifdef VERBOSE_DEBUG
#  define DEBUG(x, y) std::cout << (x) << (y) << std::endl;
#else
#  define DEBUG(x, y)
#endif

MugastIdentification::MugastIdentification() :
        lightParticles({"m1_z1", "m2_z1", "m3_z1", "m3_z2", "m4_z2"}),
        cutsMg({1, 3, 4, 5, 7, 11}),
        strips({"X", "Y"}),
        layers({"ice_front","havar_front","he3_front","he3_back","havar_back","ice_back","al_front"}),
        useConstantThickness(false),
        with_cuts(true),
        havarThickness(3.8E-3 * UNITS::mm), //in mm
        data(nullptr){}

MugastIdentification::~MugastIdentification()
{
    delete gasThickness;
    delete havarAngle;
    delete data;

    for (const auto &MG : cutsMg)
    {
        delete calibrationsTy[MG];
    }
    for (const auto &particle : lightParticles)
    {
        if (particle == "m4_z2")
            continue; //TODO: generate table
        for (const auto &layer : layers)
        {
            delete energyLoss[particle][layer];
        }
    }
    for (const auto &layer : layers)
    {
        delete energyLoss["beam"][layer];
    }
}

bool MugastIdentification::Initialize(const double &beam_energy,
                                      const TVector3 &target_pos){
    DEBUG("------------>MugastIdentification::initialize()", "");
    this->beamEnergy = beam_energy;
    this->targetPos = target_pos;

    if (useConstantThickness){
        currentIceThickness.first = 3E-2 * UNITS::mm;
        currentIceThickness.second = 3E-2 * UNITS::mm;
    }

    //Setup of reactions////////////////////////////////
    //Reactions where both fragments have been identified///////////////
    reaction.emplace("M48_Z19_m1_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(1, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(48, 19, 0, 0, 0)})));
    reaction.emplace("M47_Z19_m2_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(2, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(47, 19, 0,0, 0)})));
    reaction.emplace("M45_Z18_m4_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(4, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(45, 18, 0,0, 0)})));
    reaction.emplace("M46_Z18_m3_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0,0, 0)})));

    //Reactions where only one fragment has been identified///////////////
    reaction.emplace("m1_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(1, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(48, 19, 0, 0, 0)})));
    reaction.emplace("m2_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(2, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(47, 19, 0,0, 0)})));
    reaction.emplace("m4_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(4, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(45, 18, 0,0, 0)})));
    reaction.emplace("m3_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0,0, 0)})));

    beamRef = &reaction.at("M47_Z19_m2_z1")->GetReactionFragment(1);

    //Setting up interpolations for target deformation//////
    initializeInterpolations();

    //Cuts Initialization///////////////////////////////////
    initializeCuts();

    //Calibration initialization////////////////////////////
    initializeCalibration();

    //Energy Loss///////////////////////////////////////////
    initializeELoss();

    return true;
}

bool MugastIdentification::initializeCuts(){

    //E-TOF cuts
    std::unordered_map<std::string, TCutG *> tmp;
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
            tmp[cut_name] = cuts.at(cut_name);
        }catch (const std::out_of_range &err){
            std::cerr << "Unable to find one of the E_TOF cuts : " << cut_name << "\n";
        }
    }
    cut_type["E_TOF"] = tmp;

    //DE-E cuts
    tmp.clear();
    tmp_cut_names.clear();

    tmp_cut_names.emplace_back("DE_E_m1_z1_MM");
    tmp_cut_names.emplace_back("DE_E_m2_z1_MM");
    tmp_cut_names.emplace_back("DE_E_m3_z1_MM");
    tmp_cut_names.emplace_back("DE_E_m3_z2_MM");
    tmp_cut_names.emplace_back("DE_E_m4_z2_MM");
    for (const auto &cut_name : tmp_cut_names){
        try{
            tmp[cut_name] = cuts.at(cut_name);
        }catch (const std::out_of_range &err){
            std::cerr << "Unable to find one of the E_TOF cuts : " << cut_name << "\n";
        }
    }
    cut_type["DE_E"] = tmp;
    return true;
}

bool MugastIdentification::initializeCalibration(){
    for (const auto &MG : cutsMg){
        std::string file_name = "./Configs/Calibrations/TY_MG" + to_string(MG) + ".cal";
        ifstream file(file_name);
        if (!file){
            calibrationsTy[MG] = nullptr;
            continue;
        }
        std::cout << "Calibration :" + file_name + " found\n";
        calibrationsTy[MG] = new Calibration(file_name, 1, nStrips);
    }
    return true;
}

bool MugastIdentification::initializeELoss(){

    std::array<std::string, 7> fragments ={"m1_z1",
                                           "m2_z1",
                                           "m3_z1",
                                           "m3_z2",
                                           "m4_z2",
                                           "m46_z18",
                                           "m47_z19"};

    std::unordered_map<std::string, std::string> layer_names;
    layer_names["ice_front"]    = "Ice";
    layer_names["havar_front"]  = "Havar";
    layer_names["he3_front"]    = "Helium";
    layer_names["he3_back"]     = "Helium";
    layer_names["havar_back"]   = "Havar";
    layer_names["ice_back"]     = "Ice";
    layer_names["al_front"]     = "Al";

    std::unordered_map<std::string, std::string> fragment_names;
    fragment_names["m1_z1"]     = "H1";
    fragment_names["m2_z1"]     = "H2";
    fragment_names["m3_z1"]     = "H3";
    fragment_names["m3_z2"]     = "He3";
    fragment_names["m4_z2"]     = "He4";
    fragment_names["m46_z18"]   = "Ar46";
    fragment_names["m47_z19"]   = "K47";

    for (const auto &it_frag : fragments){
        for (const auto &it_layer : layers){
            energyLoss[it_frag][it_layer] = nullptr;
        }
    }

    std::string e_loss_file_path("./Configs/ELossTables/");

    int beamPrecision = 100;       //Beam can have lower precision
    int fragmentPrecision = 200;   //Beam can have lower precision

    if (e_loss_file_path != nullptr){
        for (const auto &it : fragments){
            int precision;
            if (it == "m46_z18")
                precision = beamPrecision;
            else
                precision = fragmentPrecision;

            for (const auto &layer : layers){
                std::string tmpString;
                tmpString = e_loss_file_path + fragment_names[it] + "_in_" + layer_names[layer] + ".txt";
                energyLoss[it][layer] = new EnergyLoss(tmpString, "SRIM", precision);
            }
        }
    }

    for (const auto &layer : layers){
        energyLoss["beam"][layer] = energyLoss["m46_z18"][layer];
    }

    currentIceThickness.first = 10E-3 * UNITS::mm; //in mm
    currentIceThickness.second = currentIceThickness.first * icePercentageSecond;
    beamEnergyMatchThreashold = 0.3 * UNITS::MeV; //In MeV

    TW_vs_ice.reserve(100);
    return true;
}

bool MugastIdentification::initializeInterpolations() {
    gasThickness = new Interpolation("./Configs/Interpolations/He_thickness_700mbar.txt");

    gasThicknessCartesian = new Interpolation("./Configs/Interpolations/He_thickness_700mbar_euc.txt");
    gasThicknessCartesian->SetUnits(UNITS::mm, UNITS::mm);

    havarAngle = new Interpolation("./Configs/Interpolations/Entrance_angle_700mbar.txt");

    auto* tmpFile = new TFile("./Configs/Interpolations/TW_Brho_M46_Z18.root");
    if (tmpFile->IsOpen()){
        TW_Brho_M46_Z18 = new Interpolation(tmpFile);
        tmpFile->Close();
        delete tmpFile;
    }else{
        std::cerr << "Unable to find TW_Brho file\n";
        TW_Brho_M46_Z18 = nullptr;
    }
    return true;
}

std::vector<std::pair<double, double>> MugastIdentification::getTWvsIce(){
    return TW_vs_ice;
}

bool MugastIdentification::Identify() {
    DEBUG("------------>MugastIdentification::identify()", "");
    //Initialization of basic structure

    for (unsigned int ii = 0; ii < fragmentCats.multiplicity; ++ii) {
        fragmentCats.Pos[ii] = TVector3((**(data->Cats)).PositionX[ii],
                                        (**(data->Cats)).PositionY[ii],
                                        (**(data->Cats)).PositionZ[ii]);
        if (fragmentCats.Pos[ii].X() == -1000) fragmentCats.Pos[ii].SetX(0);
        if (fragmentCats.Pos[ii].Y() == -1000) fragmentCats.Pos[ii].SetY(0);

        fragmentCats.Charge[ii] = TVector3((**(data->Cats)).ChargeX[ii],
                                           (**(data->Cats)).ChargeY[ii],
                                            0);

        if (ii==0){
            for (unsigned int jj=0; jj<focusScale.size(); ++jj) {
                fragmentCats.PosOnTarget[jj] = (fragmentCats.Pos[0]*focusScale[jj]);
                fragmentCats.PosOnTarget[jj].SetZ(0);
                fragmentCats.Dir[jj] = fragmentCats.Pos[0]-fragmentCats.PosOnTarget[jj];
            }
        }
    }

    for (unsigned int ii = 0; ii < fragment.multiplicity; ++ii) {
        fragment.Indentified[ii] = false;
        fragment.Pos[ii] = TVector3((**(data->Mugast)).PosX[ii],
                                    (**(data->Mugast)).PosY[ii],
                                    (**(data->Mugast)).PosZ[ii]);
        fragment.EmissionDirection[ii] = fragment.Pos[ii] - targetPos;

        if (ii==0){
            for(unsigned int jj=0; jj<beamPositions.size();++jj){
                fragment.BeamPosition_uncentered[jj] = TVector3(beamPositions[jj].first,beamPositions[jj].second, 0);
                fragment.EmissionDirection_uncentered[jj] = fragment.EmissionDirection[0] - fragment.BeamPosition_uncentered[jj];
            }
            for (unsigned int jj=0; jj<focusScale.size(); ++jj) {
                fragment.BeamPosition_corrected[jj] = fragmentCats.PosOnTarget[jj];
                fragment.EmissionDirection_corrected[jj] = fragment.EmissionDirection[0] -fragmentCats.PosOnTarget[jj];
            }
        }

        fragment.TelescopeNormal[ii] = TVector3((**(data->Mugast)).TelescopeNormalX[ii],
                                                (**(data->Mugast)).TelescopeNormalY[ii],
                                                (**(data->Mugast)).TelescopeNormalZ[ii]);
        fragment.SI_E[ii] = (**(data->Mugast)).DSSD_E[ii];
        fragment.SI_E2[ii] = (**(data->Mugast)).SecondLayer_E[ii];
        fragment.Tot_E[ii] = fragment.SI_E[ii];
        fragment.SI_X[ii] = (**(data->Mugast)).DSSD_X[ii];
        fragment.SI_Y[ii] = (**(data->Mugast)).DSSD_Y[ii];
        fragment.SI_T[ii] = (**(data->Mugast)).DSSD_T[ii];
        fragment.T2[ii] = (**(data->Mugast)).SecondLayer_T[ii];
        fragment.MG[ii] = (**(data->Mugast)).TelescopeNumber[ii];
    }
    for (unsigned int ii = 0; ii < fragmentMust.multiplicity; ++ii) {
        fragmentMust.Indentified[ii] = false;
        fragmentMust.Pos[ii] = TVector3((**(data->Must2)).PosX[ii],
                                        (**(data->Must2)).PosY[ii],
                                        (**(data->Must2)).PosZ[ii]);
        fragmentMust.EmissionDirection[ii] = fragmentMust.Pos[ii] - targetPos;


        if (ii==0){
            for(unsigned int jj=0; jj<beamPositions.size();++jj){
                fragmentMust.BeamPosition_uncentered[jj] = TVector3(beamPositions[jj].first,beamPositions[jj].second, 0);
                fragmentMust.EmissionDirection_uncentered[jj] = fragmentMust.EmissionDirection[0] - fragmentMust.BeamPosition_uncentered[jj];
            }
            for (unsigned int jj=0; jj<focusScale.size(); ++jj) {
                fragmentMust.BeamPosition_corrected[jj] = fragmentCats.PosOnTarget[jj];
                fragmentMust.EmissionDirection_corrected[jj] = fragmentMust.EmissionDirection[0] -fragmentCats.PosOnTarget[jj];
            }
        }

        fragmentMust.TelescopeNormal[ii] = TVector3(0, 0, 0);
        fragmentMust.SI_E[ii] = (**(data->Must2)).Si_E[ii];
        fragmentMust.CsI_E[ii] = (**(data->Must2)).CsI_E[ii];
        fragmentMust.CsI_T[ii] = (**(data->Must2)).CsI_T[ii];
        fragmentMust.SI_E2[ii] = (**(data->Must2)).Si_EY[ii];
        fragmentMust.Tot_E[ii] = fragmentMust.SI_E[ii] + fragmentMust.CsI_E[ii];
        fragmentMust.SI_X[ii] = (**(data->Must2)).Si_X[ii];
        fragmentMust.SI_Y[ii] = (**(data->Must2)).Si_Y[ii];
        fragmentMust.SI_T[ii] = (**(data->Must2)).Si_T[ii];
        fragmentMust.T2[ii] = (**(data->Must2)).Si_TY[ii];//TODO: check if second layers are different
        fragmentMust.MG[ii] = (**(data->Must2)).TelescopeNumber[ii];
    }

    DEBUG("------------>finished: setting up fragment", "");

    brho = TW_Brho_M46_Z18->Evaluate(**(data->TW));
    finalBeamEnergy = beamRef->GetEkFromBrho_Q(brho, chargeStateInterpolation);
    initialBeamEnergy = InitialBeamEnergy(finalBeamEnergy, currentIceThickness.first);

    if (!useConstantThickness) {
        identifyIceThickness();
    }

    for (const auto &it : reaction) {
        it.second->SetBeamEnergy(MiddleTargetBeamEnergy(finalBeamEnergy));
    }

    DEBUG("------------>starting: recalibrations", "");
    //Applying (time) re-calibrations
    for (unsigned int ii = 0; ii < fragment.multiplicity; ++ii) {
        if (calibrationsTy[fragment.MG[ii]] == nullptr)
            fragment.T[ii] = fragment.SI_T[ii];
        else
            fragment.T[ii] = calibrationsTy[fragment.MG[ii]]
                    ->Evaluate(fragment.SI_T[ii], fragment.SI_Y[ii]);
    };

    DEBUG("------------>finished: calibrations\n", "");

    //Identification with E TOF
    identifyMugast();
    reconstructEnergyMugast();
    identifyMust2();
    reconstructEnergyMust2();
    return true;
}

bool MugastIdentification::identifyMugast() {
    std::unordered_map<std::string, TCutG*>::const_iterator foundCut;

    for (unsigned int ii = 0; ii < fragment.multiplicity; ++ii) {
        fragment.Particle[ii] = "NONE";
        for (const auto &cut_it : lightParticles) {
            foundCut = cut_type["E_TOF"].find("E_TOF_" + cut_it + "_MG" + std::to_string(static_cast<int>(fragment.MG[ii])));

            if (foundCut == cut_type["E_TOF"].end()) continue;

            if (foundCut->second->IsInside(fragment.Tot_E[ii], fragment.T[ii])) {
                if (fragment.Indentified[ii])
                    throw std::runtime_error("Overlapping MUGAST E TOF gates :" +
                                             cut_it +
                                             "\tMG" +
                                             fragment.MG[ii] +
                                             "\n");

                fragment.M[ii] = std::stoi(cut_it.substr(cut_it.find_first_of('m') + 1,
                                                         cut_it.find_first_of('_') - cut_it.find_first_of('m') -
                                                         1));

                fragment.Z[ii] = std::stoi(cut_it.substr(cut_it.find_first_of('z') + 1,
                                                         cut_it.find_first_of('_') - cut_it.find_first_of('z') -
                                                         1));

                fragment.Indentified[ii] = true;
                fragment.Particle[ii] = cut_it;
            }
        }
        if (!fragment.Indentified[ii]) {
            fragment.M[ii] = 0;
            fragment.Z[ii] = 0;
        }
    };
    DEBUG("------------>finished: searching cuts", "");
    return true;
}

bool MugastIdentification::identifyMust2() {
    std::unordered_map<std::string, TCutG*>::const_iterator foundCut;

    for (unsigned int ii = 0; ii < fragmentMust.multiplicity; ++ii) {
        fragmentMust.Particle[ii] = "NONE";
        for (const auto &cut_it : lightParticles) {
            foundCut = cut_type["DE_E"].find("DE_E_" + cut_it + "_MM");

            if (foundCut == cut_type["DE_E"].end()) continue;

            if (foundCut->second->IsInside(fragmentMust.Tot_E[ii], fragmentMust.SI_E[ii])) {
                if (fragmentMust.Indentified[ii]) continue;
//TODO:                    throw std::runtime_error("Overlapping MUST2 DE E gates :" +
//TODO:                                             cut_it +
//TODO:                                             "\n");


                fragmentMust.M[ii] = std::stoi(cut_it.substr(cut_it.find_first_of('m') + 1,
                                                             cut_it.find_first_of('_') - cut_it.find_first_of('m') -
                                                             1));

                fragmentMust.Z[ii] = std::stoi(cut_it.substr(cut_it.find_first_of('z') + 1,
                                                             cut_it.find_first_of('_') - cut_it.find_first_of('z') -
                                                             1));

                fragmentMust.Indentified[ii] = true;
                fragmentMust.Particle[ii] = cut_it;
            }
        }
        if (!fragmentMust.Indentified[ii]) {
            fragmentMust.M[ii] = 0;
            fragmentMust.Z[ii] = 0;
        }
    };
    DEBUG("------------>finished: searching cuts", "");
    return true;
}

bool MugastIdentification::reconstructEnergyMugast() {
    //Energy reconstruction
    std::unordered_map<std::string,std::unordered_map<std::string, EnergyLoss *>>::iterator eloss_it;
    for (unsigned int ii = 0; ii < fragment.multiplicity; ++ii){

        //If not identified, no need to compute anything
        if (!fragment.Indentified[ii]){
            fragment.E[ii] = fragment.Tot_E[ii];
            fragment.Ex[ii] = -1000;
            continue;
        }

        eloss_it = energyLoss.find("m" +
                                   std::to_string(fragment.M[ii]) +
                                   "_z" +
                                   std::to_string(fragment.Z[ii]));
        if (eloss_it == energyLoss.end())
            return false;

        double theta{fragment.EmissionDirection[ii].Angle(TVector3(0, 0, -1))};

        //Passivation layer
        double tmpEn{0};
        double tmpEn2{0};

        tmpEn = eloss_it->second["al_front"]->EvaluateInitialEnergy(fragment.Tot_E[ii],
                                                                    0.4E-3*UNITS::mm, //Units in mm!
                                                                    fragment.EmissionDirection[ii]
                                                                            .Angle(fragment.TelescopeNormal[ii]));//TODO: telescope normal probabily not correct

        tmpEn = eloss_it->second["ice_front"]->EvaluateInitialEnergy(tmpEn,
                                                                     currentIceThickness.first,
                                                                     havarAngle->Evaluate(theta));

        tmpEn = eloss_it->second["havar_front"]->EvaluateInitialEnergy(tmpEn,
                                                                       havarThickness,
                                                                       havarAngle->Evaluate(theta));

        tmpEn2 = tmpEn;

        tmpEn = eloss_it->second["he3_front"]->EvaluateInitialEnergy(tmpEn,
                                                                     gasThickness->Evaluate(theta) * UNITS::mm,
                                                                     0.);

        fragment.E[ii] = tmpEn;

        if (1 || ii==0) {
            for (unsigned int jj = 0; jj < beamPositions.size(); ++jj) {
                try {
                    fragment.E_uncentered[ii].push_back(eloss_it->second["he3_front"]->EvaluateInitialEnergy(tmpEn2,
                                                                                                             computeDistanceInGas(
                                                                                                                     fragment.EmissionDirection_uncentered[jj],
                                                                                                                     fragment.BeamPosition_uncentered[jj]) *UNITS::mm,
                                                                                                             0)
                    );
                } catch (std::runtime_error &e) {
                    std::cerr << "Unable to find correct un-centered distance (Mugast) : " << e.what() << std::endl;
                    fragment.E_uncentered[ii].push_back(tmpEn);
                }
            }
            for (unsigned int jj = 0; jj < focusScale.size(); ++jj) {
                try {
                    fragment.E_corrected[ii].push_back(eloss_it->second["he3_front"]->EvaluateInitialEnergy(tmpEn2,
                                                                                                             computeDistanceInGas(
                                                                                                                     fragment.EmissionDirection_corrected[jj],
                                                                                                                     fragmentCats.PosOnTarget[jj]) *UNITS::mm,
                                                                                                             0)
                    );
                } catch (std::runtime_error &e) {
                    std::cerr << "Unable to find correct un-centered distance (Mugast) : " << e.what() << std::endl;
                    fragment.E_uncentered[ii].push_back(tmpEn);
                }
            }
        }

        if ((reactionIt = reaction.find("M" + std::to_string(data->VAMOS_id_M) +
                                        "_Z" + std::to_string(data->VAMOS_id_Z) +
                                        "_" +
                                        fragment.Particle[ii])) != reaction.end()){
            //Vamos and mugast fragment found
            fragment.Ex[ii] = reactionIt->second->Set_Ek_Theta(fragment.E[ii], fragment.EmissionDirection[ii].Theta());//WARNING this could be wring
            fragment.Ex_uncentered[ii].reserve(beamPositions.size());
            fragment.Ex_corrected[ii].reserve(focusScale.size());
            if (1 || ii == 0) {
                for (unsigned int jj = 0; jj < beamPositions.size(); ++jj) {
                    fragment.Ex_uncentered[ii].push_back(reactionIt->second->Set_Ek_Theta(fragment.E_uncentered[ii][jj],
                                                                                          fragment.EmissionDirection_uncentered[jj].Theta()));
                }
                for (unsigned int jj = 0; jj < focusScale.size(); ++jj) {
                    fragment.Ex_corrected[ii].push_back(reactionIt->second->Set_Ek_Theta(fragment.E_corrected[ii][jj],
                                                                                          fragment.EmissionDirection_corrected[jj].Theta()));
                }
            }

            fragment.E_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_Ek_cm();
        }
        else{
            if ((reactionIt = reaction.find(fragment.Particle[ii])) != reaction.end()) {
                //Only mugast fragment found
                fragment.Ex[ii] = reactionIt->second->Set_Ek_Theta(fragment.E[ii],
                                                                   fragment.EmissionDirection[ii].Theta());//WARNING this could be wring
                fragment.E_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_Ek_cm();

            }else
                fragment.Ex[ii] = -2000;
        }
    }
    DEBUG("------------>Finished : MugastIdentification::identify()", "");
    return true;
}

bool MugastIdentification::reconstructEnergyMust2() {
    //Energy reconstruction
    std::unordered_map<std::string, std::unordered_map<std::string, EnergyLoss *>>::iterator eloss_it;
    for (unsigned int ii = 0; ii < fragmentMust.multiplicity; ++ii) {

        //If not identified, no need to compute anything
        if (!fragmentMust.Indentified[ii]) {
            fragmentMust.E[ii] = fragmentMust.Tot_E[ii];
            fragmentMust.Ex[ii] = -1000;
            continue;
        }

        eloss_it = energyLoss.find("m" +
                                   std::to_string(fragmentMust.M[ii]) +
                                   "_z" +
                                   std::to_string(fragmentMust.Z[ii]));
        if (eloss_it == energyLoss.end())
            return false;

        double theta{fragmentMust.EmissionDirection[ii].Angle(TVector3(0, 0, 1))};

        //Passivation layer
        double tmpEn{0};
        double tmpEn2{0};

        tmpEn = eloss_it->second["al_front"]->EvaluateInitialEnergy(fragmentMust.Tot_E[ii],
                                                                    0.4E-3 * UNITS::mm, //Units in mm!
                                                                    0);//TODO: add telescope normal

        tmpEn = eloss_it->second["ice_front"]->EvaluateInitialEnergy(tmpEn,
                                                                     currentIceThickness.first * icePercentageSecond,
                                                                     havarAngle->Evaluate(theta));

        tmpEn = eloss_it->second["havar_front"]->EvaluateInitialEnergy(tmpEn,
                                                                       havarThickness,
                                                                       havarAngle->Evaluate(theta));

        tmpEn2 = tmpEn;

        tmpEn = eloss_it->second["he3_front"]->EvaluateInitialEnergy(tmpEn,
                                                                     gasThickness->Evaluate(theta) * UNITS::mm,
                                                                     0.);

        fragmentMust.E[ii] = tmpEn;

        if (1 || ii==0) {
            for (unsigned int jj = 0; jj < beamPositions.size(); ++jj) {
                try {
                    fragmentMust.E_uncentered[ii].push_back(eloss_it->second["he3_front"]->EvaluateInitialEnergy(tmpEn2,
                                                                                                             computeDistanceInGas(
                                                                                                                     fragmentMust.EmissionDirection_uncentered[jj],
                                                                                                                     fragmentMust.BeamPosition_uncentered[jj]) *UNITS::mm,
                                                                                                             0)
                    );
                } catch (std::runtime_error &e) {
                    std::cerr << "Unable to find correct un-centered distance (Mugast) : " << e.what() << std::endl;
                    fragmentMust.E_uncentered[ii].push_back(tmpEn);
                }
            }
            for (unsigned int jj = 0; jj < focusScale.size(); ++jj) {
                try {
                            fragmentMust.E_corrected[ii].push_back(eloss_it->second["he3_front"]->EvaluateInitialEnergy(tmpEn2,
                                                                                                                        computeDistanceInGas(
                                                                                                                    fragmentMust.EmissionDirection_corrected[jj],
                                                                                                                    fragmentCats.PosOnTarget[jj]) *UNITS::mm,
                                                                                                            0)
                    );
                } catch (std::runtime_error &e) {
                    std::cerr << "Unable to find correct un-centered distance (Mugast) : " << e.what() << std::endl;
                    fragmentMust.E_uncentered[ii].push_back(tmpEn);
                }
            }
        }

        if ((reactionIt = reaction.find("M" + std::to_string(data->VAMOS_id_M) +
                                        "_Z" + std::to_string(data->VAMOS_id_Z) +
                                        "_" +
                                        fragmentMust.Particle[ii])) != reaction.end()) {
            //Vamos and Must2 fragments found
            fragmentMust.Ex[ii] = reactionIt->second->Set_Ek_Theta(fragmentMust.E[ii],
                                                                   fragmentMust.EmissionDirection[ii].Theta());

            fragmentMust.Ex_uncentered[ii].reserve(beamPositions.size());
            fragmentMust.Ex_corrected[ii].reserve(focusScale.size());
            if (1 || ii == 0) {
                for (unsigned int jj = 0; jj < beamPositions.size(); ++jj) {
                    fragmentMust.Ex_uncentered[ii].push_back(reactionIt->second->Set_Ek_Theta(fragmentMust.E_uncentered[ii][jj],
                                                                                          fragmentMust.EmissionDirection_uncentered[jj].Theta()));
                }
                for (unsigned int jj = 0; jj < focusScale.size(); ++jj) {
                    fragmentMust.Ex_corrected[ii].push_back(reactionIt->second->Set_Ek_Theta(fragmentMust.E_corrected[ii][jj],
                                                                                          fragmentMust.EmissionDirection_corrected[jj].Theta()));
                }
            }

            fragmentMust.E_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_Ek_cm();
        } else {
            if ((reactionIt = reaction.find(fragmentMust.Particle[ii])) != reaction.end()) {
                //Only Must2 fragment found
                fragmentMust.Ex[ii] = reactionIt->second->Set_Ek_Theta(fragmentMust.E[ii],
                                                                       fragmentMust.EmissionDirection[ii].Theta());
                fragmentMust.E_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_Ek_cm();
            }else
                fragmentMust.Ex[ii] = -2000;
        }
    }
    DEBUG("------------>Finished : MugastIdentification::identify()", "");
    return true;
}

void MugastIdentification::identifyIceThickness(){
    DEBUG("------------>Started: MugastIdentification::identifyIceThickness()", "");
    DEBUG("Final beam energy : " , finalBeamEnergy);
    DEBUG("Computed initial beam energy : " , initialBeamEnergy);
    DEBUG(" With thickness : " , currentIceThickness.first);
    DEBUG(" With brho : " , brho);
    DEBUG(" and mass : " , beamRef->Get_M2());

    if (abs(beamEnergy - initialBeamEnergy) > beamEnergyMatchThreashold){
        DEBUG("Before minimizer call, ice_thickness.first = ", currentIceThickness.first);
        double tmp_threashold{2E-4*UNITS::mm};
        double tmp_precision{0.1*UNITS::MeV};

        iceThicknessMinimizer.reset(new Minimizer([this](const double &tck) {
                                                      return pow(this->beamEnergy
                                                                 - this->InitialBeamEnergy(this->finalBeamEnergy, tck), 2)
                                                             /(this->beamEnergy * this->beamEnergy);
                                                  },
                                                  currentIceThickness.first,        //starting value
                                                  2E-4 * currentIceThickness.first, //learning rate
                                                  tmp_threashold,                     //threashold
                                                  100,                                //max_steps
                                                  1,                                  //quenching
                                                  1.E-4 * currentIceThickness.first)  //h
        );

        while (fabs(beamEnergy - InitialBeamEnergy(finalBeamEnergy, currentIceThickness.first)) > tmp_precision){
            tmp_threashold/=2;
            iceThicknessMinimizer->SetThreashold(tmp_threashold);
            currentIceThickness.first = iceThicknessMinimizer->Minimize();
            currentIceThickness.second = currentIceThickness.first * icePercentageSecond;
        }

        DEBUG("After minimizer call, ice_thickness.first = ", currentIceThickness.first);
        DEBUG(" energy difference : ", beamEnergy - InitialBeamEnergy(finalBeamEnergy, currentIceThickness.first));
    }
    DEBUG("------------>Finished: MugastIdentification::identifyIceThickness()", "");
}

double MugastIdentification::InitialBeamEnergy(double beam_energy_from_brho, double tck){
    DEBUG("------------>Started: MugastIdentification::InitialBeamEnergy(), thickness : ", tck);
    beam_energy_from_brho =
            energyLoss["beam"]["ice_back"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            tck * icePercentageSecond,
                                            0);

    beam_energy_from_brho =
            energyLoss["beam"]["havar_back"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            havarThickness,
                                            0);

    beam_energy_from_brho =
            energyLoss["beam"]["he3_front"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            averageBeamThickness,
                                            0);

    beam_energy_from_brho =
            energyLoss["beam"]["havar_front"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            havarThickness,
                                            0);

    beam_energy_from_brho =
            energyLoss["beam"]["ice_front"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            tck,
                                            0);
    DEBUG("------------>Finished: MugastIdentification::InitialBeamEnergy()", "");
    return beam_energy_from_brho;
}

double MugastIdentification::MiddleTargetBeamEnergy(double beam_energy_from_brho){
    beam_energy_from_brho =
            energyLoss["beam"]["ice_back"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            currentIceThickness.first,
                                            0);

    beam_energy_from_brho =
            energyLoss["beam"]["havar_back"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            havarThickness,
                                            0);

    beam_energy_from_brho =
            energyLoss["beam"]["he3_front"]
                    ->EvaluateInitialEnergy(beam_energy_from_brho,
                                            averageBeamThickness / 2.,
                                            0);

    return beam_energy_from_brho;
}

double MugastIdentification::computeDistanceInGas(TVector3 vectDet, const TVector3& incidentPos) {

    double tmp_threashold{1E-3*UNITS::mm};

    if (vectDet.Z()>0){
        distanceMinimizer.reset(new Minimizer([&,this](double x) {
                                                  auto vectTmp = x*vectDet+incidentPos;
                                                  return pow(vectTmp.Z()-this->gasThicknessCartesian->Evaluate(sqrt(pow(vectTmp.X(), 2) + pow(vectTmp.Y(), 2))),
                                                             2);
                                              },
                                              1.5*UNITS::mm/vectDet.CosTheta(),   //starting value
                                              1E-2*UNITS::mm,                     //learning rate
                                              tmp_threashold,                     //threashold
                                              100,                                //max_steps
                                              1,                                  //quenching
                                              1E-2*UNITS::mm)                     //h
        );
    }else{
        distanceMinimizer.reset(new Minimizer([&,this](double x) {
                                                  //auto vectTmp = x*vectDet+incidentPos;
                                                  auto vectTmp = x*vectDet+incidentPos;
                                                  return pow(vectTmp.Z()+this->gasThicknessCartesian->Evaluate(sqrt(pow(vectTmp.X(), 2) + pow(vectTmp.Y(), 2))),
                                                             2);
                                              },
                                              1.5*UNITS::mm/vectDet.CosTheta(),   //starting value
                                              1E-2*UNITS::mm,                     //learning rate
                                              tmp_threashold,                     //threashold
                                              100,                                //max_steps
                                              1,                                  //quenching
                                              1E-2*UNITS::mm)                     //h
        );
    }

    distanceMinimizer->Minimize();
    auto distance = ((distanceMinimizer->GetValue())*vectDet + incidentPos).Mag();

    if (distance > 100*UNITS::m) throw std::runtime_error("Distance too high\n");

    return distance;
}