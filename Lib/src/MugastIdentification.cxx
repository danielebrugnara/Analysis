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
        havarThickness(3.8 * UNITS::micrometer), //in mm
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

bool MugastIdentification::initialize(const double &beamEnergy,
                                      const TVector3 &targetPos){
    DEBUG("------------>MugastIdentification::initialize()", "");
    this->beamEnergy = beamEnergy;
    this->targetPos = targetPos;

    if (useConstantThickness){
        currentIceThickness.first = 3E-2 * UNITS::mm;
        currentIceThickness.second = 3E-2 * UNITS::mm;
    }

    //Setup of reactions////////////////////////////////
    //Reactions where both fragments have been identified///////////////
    reaction.emplace("M48_Z19_m1_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(1, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(48, 19, 0, 0, 0)})));
    reaction.emplace("M47_Z19_m2_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(2, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(47, 19, 0,0, 0)})));
    reaction.emplace("M45_Z18_m4_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(4, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(45, 18, 0,0, 0)})));
    reaction.emplace("M46_Z18_m3_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0,0, 0)})));

    //Reactions where only one fragment has been identified///////////////
    reaction.emplace("m1_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(1, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(48, 19, 0, 0, 0)})));
    reaction.emplace("m2_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(2, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(47, 19, 0,0, 0)})));
    reaction.emplace("m4_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(4, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(45, 18, 0,0, 0)})));
    reaction.emplace("m3_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beamEnergy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0,0, 0)})));

    for (const auto& it: reaction){
        it.second->ChooseFixed(3);
        //it.second->ChooseExFixed(4);
        it.second->ChooseExFixed(3);
    }

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
    std::vector<std::string> tmpCutNames;

    tmpCutNames.emplace_back("E_TOF_m1_z1_MG1");
    tmpCutNames.emplace_back("E_TOF_m1_z1_MG3");
    tmpCutNames.emplace_back("E_TOF_m1_z1_MG4");
    tmpCutNames.emplace_back("E_TOF_m1_z1_MG5");
    tmpCutNames.emplace_back("E_TOF_m1_z1_MG7");
    tmpCutNames.emplace_back("E_TOF_m1_z1_MG11");

    tmpCutNames.emplace_back("E_TOF_m2_z1_MG1");
    tmpCutNames.emplace_back("E_TOF_m2_z1_MG3");
    tmpCutNames.emplace_back("E_TOF_m2_z1_MG4");
    tmpCutNames.emplace_back("E_TOF_m2_z1_MG5");
    tmpCutNames.emplace_back("E_TOF_m2_z1_MG7");
    tmpCutNames.emplace_back("E_TOF_m2_z1_MG11");

    tmpCutNames.emplace_back("E_TOF_m4_z2_MG1");
    tmpCutNames.emplace_back("E_TOF_m4_z2_MG3");
    tmpCutNames.emplace_back("E_TOF_m4_z2_MG4");
    tmpCutNames.emplace_back("E_TOF_m4_z2_MG5");
    tmpCutNames.emplace_back("E_TOF_m4_z2_MG7");
    tmpCutNames.emplace_back("E_TOF_m4_z2_MG11");

    for (const auto &cut_name : tmpCutNames){
        try{
            tmp[cut_name] = cuts.at(cut_name);
        }catch (const std::out_of_range &err){
            std::cerr << "Unable to find one of the E_TOF cuts : " << cut_name << "\n";
        }
    }
    cut_type["E_TOF"] = tmp;

    //DE-E cuts
    tmp.clear();
    tmpCutNames.clear();

    tmpCutNames.emplace_back("DE_E_m1_z1_MM");
    tmpCutNames.emplace_back("DE_E_m2_z1_MM");
    tmpCutNames.emplace_back("DE_E_m3_z1_MM");
    tmpCutNames.emplace_back("DE_E_m3_z2_MM");
    tmpCutNames.emplace_back("DE_E_m4_z2_MM");
    for (const auto &cut_name : tmpCutNames){
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

    std::unordered_map<std::string, std::string> layerNames;
    layerNames["ice_front"]    = "Ice";
    layerNames["havar_front"]  = "Havar";
    layerNames["he3_front"]    = "Helium";
    layerNames["he3_back"]     = "Helium";
    layerNames["havar_back"]   = "Havar";
    layerNames["ice_back"]     = "Ice";
    layerNames["al_front"]     = "Al";

    std::unordered_map<std::string, std::string> fragmentNames;
    fragmentNames["m1_z1"]     = "H1";
    fragmentNames["m2_z1"]     = "H2";
    fragmentNames["m3_z1"]     = "H3";
    fragmentNames["m3_z2"]     = "He3";
    fragmentNames["m4_z2"]     = "He4";
    fragmentNames["m46_z18"]   = "Ar46";
    fragmentNames["m47_z19"]   = "K47";

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
                tmpString = e_loss_file_path + fragmentNames[it] + "_in_" + layerNames[layer] + ".txt";
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

    //TW_vs_ice.reserve(100);
    return true;
}

bool MugastIdentification::initializeInterpolations() {
    gasThickness = new Interpolation("./Configs/Interpolations/He_thickness_700mbar.txt");

    gasThicknessCartesian = new Interpolation("./Configs/Interpolations/He_thickness_700mbar_euc.txt");
    gasThicknessDerivativeCartesian = new Interpolation("./Configs/Interpolations/He_thickness_derivative_700mbar_euc.txt");
    //gasThicknessCartesian->SetUnits(UNITS::mm, UNITS::mm);

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

//std::vector<std::pair<double, double>> MugastIdentification::getTWvsIce(){
//    return TW_vs_ice;
//}

template<class D>
bool MugastIdentification::fillInitialData(MugastData & localFragment, const D* nplData) {

    for (unsigned int ii = 0; ii < localFragment.multiplicity; ++ii) {
        localFragment.Indentified[ii] = false;
        localFragment.Pos[ii]               = TVector3(nplData->PosX[ii]*UNITS::mm,nplData->PosY[ii]*UNITS::mm,nplData->PosZ[ii]*UNITS::mm);
        localFragment.EmissionDirection[ii] = localFragment.Pos[ii] - targetPos;

        if (ii==0){
            for(unsigned int jj=0; jj<beamPositions.size();++jj){
                localFragment.BeamPosition_uncentered[jj]       = TVector3(beamPositions[jj].first, beamPositions[jj].second, 0);
                localFragment.EmissionDirection_uncentered[jj]  = localFragment.EmissionDirection[0] - localFragment.BeamPosition_uncentered[jj];
            }
            for (unsigned int jj=0; jj<focusScale.size(); ++jj) {
                localFragment.BeamPosition_corrected[jj]        = fragmentCats.PosOnTarget[jj];
                localFragment.EmissionDirection_corrected[jj]   = localFragment.EmissionDirection[0] - fragmentCats.PosOnTarget[jj];
                localFragment.BeamDirection_corrected[jj]       = localFragment.BeamPosition_corrected[jj] - fragmentCats.Pos[0];
            }
        }

        localFragment.TelescopeNormal[ii] = TVector3(0,0,0);
        localFragment.MG[ii] = nplData->TelescopeNumber[ii];
        auto* mugast = dynamic_cast<const TMugastPhysics*>(nplData);
        auto* must2 = dynamic_cast<const TMust2Physics*>(nplData);
        if (mugast != nullptr && must2 == nullptr){
            localFragment.TelescopeNormal[ii] = TVector3(   mugast->TelescopeNormalX[ii],
                                                            mugast->TelescopeNormalY[ii],
                                                            mugast->TelescopeNormalZ[ii]);
            localFragment.SI_E[ii]  = mugast->DSSD_E[ii]*UNITS::MeV;
            localFragment.SI_E2[ii] = mugast->SecondLayer_E[ii]*UNITS::MeV;
            localFragment.Tot_E[ii] = localFragment.SI_E[ii]*UNITS::MeV;
            localFragment.SI_X[ii]  = mugast->DSSD_X[ii];
            localFragment.SI_Y[ii]  = mugast->DSSD_Y[ii];
            localFragment.SI_T[ii]  = mugast->DSSD_T[ii];
            localFragment.T2[ii]    = mugast->SecondLayer_T[ii];
        }else if (mugast == nullptr && must2 != nullptr){
            localFragment.TelescopeNormal[ii] = TVector3(0, 0, 0);
            localFragment.SI_E[ii]  = must2->Si_E[ii]*UNITS::MeV;
            localFragment.SI_E2[ii] = must2->Si_EY[ii]*UNITS::MeV;
            //Tot_E filled after calling this function
            localFragment.SI_X[ii]  = must2->Si_X[ii];
            localFragment.SI_Y[ii]  = must2->Si_Y[ii];
            localFragment.SI_T[ii]  = must2->Si_T[ii];
            localFragment.T2[ii]    = must2->Si_TY[ii];//TODO: check if second layers are different
        }else{
            throw std::runtime_error("Something wrong with pointers\n");
        }
    }
    DEBUG("------------>starting: recalibrations", "");
    //Applying (time) re-calibrations
    for (unsigned int ii = 0; ii < localFragment.multiplicity; ++ii) {
        if (calibrationsTy[localFragment.MG[ii]] == nullptr)
            localFragment.T[ii] = localFragment.SI_T[ii];
        else
            localFragment.T[ii] = calibrationsTy[localFragment.MG[ii]]
                    ->Evaluate(localFragment.SI_T[ii], localFragment.SI_Y[ii]);
    };
    DEBUG("------------>finished: calibrations\n", "");
    return true;
}

bool MugastIdentification::identify() {
    DEBUG("------------>MugastIdentification::identify()", "");
    //Initialization of basic structure
    //Cats2
    for (unsigned int ii = 0; ii < fragmentCats.multiplicity; ++ii) {
        if ((**(data->Cats)).PositionX.empty()) { //Assume beam in the center of the target
            fragmentCats.Pos[0] = TVector3(0,0,-2045*UNITS::mm);
            fragmentCats.Charge[ii] = TVector3(0,0,0);
        }else{
            fragmentCats.Pos[ii] = TVector3((**(data->Cats)).PositionX[ii]*UNITS::mm,
                                            (**(data->Cats)).PositionY[ii]*UNITS::mm,
                                            (**(data->Cats)).PositionZ[ii]*UNITS::mm);
            if (fragmentCats.Pos[ii].X() == -1000) fragmentCats.Pos[ii].SetX(0);
            if (fragmentCats.Pos[ii].Y() == -1000) fragmentCats.Pos[ii].SetY(0);
            if (fragmentCats.Pos[ii].Y() == -1000) fragmentCats.Pos[ii].SetZ(-2045*UNITS::mm);

            fragmentCats.Charge[ii] = TVector3((**(data->Cats)).ChargeX[ii],
                                               (**(data->Cats)).ChargeY[ii],
                                               0);
        }
        if (ii==0){
            for (unsigned int jj=0; jj<focusScale.size(); ++jj) {
                fragmentCats.PosOnTarget[jj] = (fragmentCats.Pos[0]*focusScale[jj]);
                fragmentCats.PosOnTarget[jj].SetZ(0);
                if (fragmentCats.PosOnTarget[jj].Mag()>targetRadius){
                    fragmentCats.PosOnTarget[jj].SetMag(targetRadius);
                }
                fragmentCats.Dir[jj] = fragmentCats.Pos[0]-fragmentCats.PosOnTarget[jj];
            }
        }
    }

    //Mugast
    fillInitialData<TMugastPhysics>(fragment, &(**(data->Mugast))); //Ugly de-referencing is ROOT's fault!

    //Must2
    fillInitialData<TMust2Physics>(fragmentMust, &(**(data->Must2))); //Ugly de-referencing is ROOT's fault!

    for (unsigned int ii = 0; ii < fragmentMust.multiplicity; ++ii) { //Things in Must2 but not in Mugast
        fragmentMust.CsI_E[ii] = (**(data->Must2)).CsI_E[ii]*UNITS::MeV;
        fragmentMust.CsI_T[ii] = (**(data->Must2)).CsI_T[ii];
        fragmentMust.Tot_E[ii] = fragmentMust.SI_E[ii] + fragmentMust.CsI_E[ii];
    }

    DEBUG("------------>finished: setting up fragment", "");

    brho = TW_Brho_M46_Z18->Evaluate(**(data->TW));
    finalBeamEnergy = beamRef->GetEkFromBrho_Q(brho, chargeStateInterpolation);
    initialBeamEnergy = InitialBeamEnergy(finalBeamEnergy, currentIceThickness.first);

    if (!useConstantThickness) {
        identifyIceThickness();
    }

    currentMidTargetBeamEnergy = middleTargetBeamEnergy(finalBeamEnergy);
    for (const auto &it : reaction) {
        it.second->SetBeamEnergy(currentMidTargetBeamEnergy);
    }

    //Identification with E TOF
    identifyMugast();
    reconstructEnergy(fragment);
    identifyMust2();
    reconstructEnergy(fragmentMust);
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

bool MugastIdentification::reconstructEnergy(MugastData & localFragment) {
    std::unordered_map<std::string,std::unordered_map<std::string, EnergyLoss *>>::iterator elossIt;
    for (unsigned int ii = 0; ii < localFragment.multiplicity; ++ii){//loop over multiplicity

        //If not identified, no need to compute anything
        if (!localFragment.Indentified[ii]){
            localFragment.E[ii] = localFragment.Tot_E[ii];
            localFragment.Ex[ii] = -1000;
            continue;
        }

        elossIt = energyLoss.find("m" +
                                  std::to_string(localFragment.M[ii]) +
                                  "_z" +
                                  std::to_string(localFragment.Z[ii]));
        if (elossIt == energyLoss.end())
            return false;

        double theta{0};

        if (localFragment.EmissionDirection[ii].Z()>0) //Must2
            theta = localFragment.EmissionDirection[ii].Angle(TVector3(0, 0, 1));
        else //Mugast
            theta = localFragment.EmissionDirection[ii].Angle(TVector3(0, 0, -1));

        //Passivation layer
        double tmpEn{0};
        double energyAfterTarget{0};


        double telescopeEntranceAngle = 0;
        if (abs(localFragment.TelescopeNormal[ii].Mag())> 0.001) //If vector is not zero
            telescopeEntranceAngle = localFragment.EmissionDirection[ii].Angle(localFragment.TelescopeNormal[ii]);

        tmpEn = elossIt->second["al_front"]->EvaluateInitialEnergy(localFragment.Tot_E[ii],
                                                                   0.4E-3*UNITS::mm, //Units in mm!
                                                                   telescopeEntranceAngle);//TODO: telescope normal probabily not correct

        energyAfterTarget = tmpEn;

        if (ii==0) {
            localFragment.EBeforeLayer[3] = tmpEn;
            localFragment.ThickessInLayer[3] = 0.4E-3 * UNITS::mm / cos(telescopeEntranceAngle);
        }

        //Deformed Target
        tmpEn = elossIt->second["ice_front"]->EvaluateInitialEnergy(energyAfterTarget,
                                                                    currentIceThickness.first,
                                                                    havarAngle->Evaluate(theta));

        if (ii==0) {
            localFragment.EBeforeLayer[2] = tmpEn;
            localFragment.ThickessInLayer[2] = currentIceThickness.first / cos(havarAngle->Evaluate(theta));
        }

        tmpEn = elossIt->second["havar_front"]->EvaluateInitialEnergy(tmpEn,
                                                                      havarThickness,
                                                                      havarAngle->Evaluate(theta));

        if (ii==0) {
            localFragment.EBeforeLayer[1] = tmpEn;
            localFragment.ThickessInLayer[1] = havarThickness / cos(havarAngle->Evaluate(theta));
        }

        tmpEn = elossIt->second["he3_front"]->EvaluateInitialEnergy(tmpEn,
                                                                    gasThickness->Evaluate(theta) * UNITS::mm,
                                                                    0.);

        if (ii==0) {
                localFragment.EBeforeLayer[0] = tmpEn;
                localFragment.ThickessInLayer[0] = gasThickness->Evaluate(theta) * UNITS::mm / cos(0);
        }

        localFragment.E[ii] = tmpEn;

        //Flat Target
        tmpEn = elossIt->second["ice_front"]->EvaluateInitialEnergy(energyAfterTarget,
                                                                    currentIceThickness.first,
                                                                    theta);

        tmpEn = elossIt->second["havar_front"]->EvaluateInitialEnergy(tmpEn,
                                                                      havarThickness,
                                                                      theta);


        tmpEn = elossIt->second["he3_front"]->EvaluateInitialEnergy(tmpEn,
                                                                    3. * UNITS::mm,
                                                                    theta);

        localFragment.EFlatTarget[ii] = tmpEn;

        //Spherical Target
        //tmpEn = elossIt->second["ice_front"]->EvaluateInitialEnergy(energyAfterTarget,
        //                                                            currentIceThickness.first,
        //                                                            theta);

        //tmpEn = elossIt->second["havar_front"]->EvaluateInitialEnergy(tmpEn,
        //                                                              havarThickness,
        //                                                              theta);


        //tmpEn = elossIt->second["he3_front"]->EvaluateInitialEnergy(tmpEn,
        //                                                            3. * UNITS::mm,
        //                                                            theta);

        localFragment.ESphericalTarget[ii] = tmpEn;

        bool checlCollision{true};
        if(checlCollision){
            TVector3 collisionVec = computeCollision(localFragment.EmissionDirection[ii],TVector3(0,0,0));
            double tentativeAngle = computeAngleOfIncidence(collisionVec, TVector3(0, 0, 0));
            if (abs(gasThickness->Evaluate(theta)*UNITS::mm-collisionVec.Mag())>0.05*UNITS::mm){
                std::cerr   << "-------------------------------\n"
                            << "should be thickness : " << gasThickness->Evaluate(theta)*UNITS::mm << std::endl
                            << "instead thickness is : " << collisionVec.Mag() << std::endl;
                //throw std::runtime_error("Something wrong in computin the collision\n");
            }
            if (abs(havarAngle->Evaluate(theta)-tentativeAngle)>1E-3){
                std::cerr   << "-------------------------------\n"
                            << "should be angle : " << havarAngle->Evaluate(theta) << std::endl
                            << "instead angle is : " << tentativeAngle << std::endl;
                //throw std::runtime_error("Something wrong in computing the angle\n");
            }
        }

        if (ii==0) {
            for (unsigned int jj = 0; jj < beamPositions.size(); ++jj) {//loop over positions
                try {
                    TVector3 collisionVec = computeCollision(localFragment.EmissionDirection_uncentered[jj],
                                                             localFragment.BeamPosition_uncentered[jj]);

                    double incidenceAngle = computeAngleOfIncidence(collisionVec,
                                                                    localFragment.BeamPosition_uncentered[jj]);

                    tmpEn = elossIt->second["ice_front"]->EvaluateInitialEnergy(energyAfterTarget,
                                                                                currentIceThickness.first,
                                                                                incidenceAngle);

                    tmpEn = elossIt->second["havar_front"]->EvaluateInitialEnergy(tmpEn,
                                                                                  havarThickness,
                                                                                  incidenceAngle);

                    localFragment.E_uncentered[ii].push_back(elossIt->second["he3_front"]->
                            EvaluateInitialEnergy(  tmpEn,
                                                    collisionVec.Mag(),
                                                    0));
                } catch (std::runtime_error &e) {
                    std::cerr << "Unable to find correct un-centered distance (Mugast) : " << e.what() << std::endl;
                    localFragment.E_uncentered[ii].push_back(tmpEn);
                    localFragment.BeamPosition_uncentered[jj] = TVector3(0,0,0);
                }
            }
            for (unsigned int jj = 0; jj < focusScale.size(); ++jj) {//loop over different focus coefficients
                try {
                    TVector3 collisionVec = computeCollision(localFragment.EmissionDirection_corrected[jj],
                                                             localFragment.BeamPosition_corrected[jj]);

                    double incidenceAngle = computeAngleOfIncidence(collisionVec,
                                                                    localFragment.BeamPosition_corrected[jj]);

                    tmpEn = elossIt->second["ice_front"]->EvaluateInitialEnergy(energyAfterTarget,
                                                                                currentIceThickness.first,
                                                                                incidenceAngle);

                    tmpEn = elossIt->second["havar_front"]->EvaluateInitialEnergy(tmpEn,
                                                                                  havarThickness,
                                                                                  incidenceAngle);

                    localFragment.E_corrected[ii].push_back(elossIt->second["he3_front"]->
                            EvaluateInitialEnergy(  tmpEn,
                                                    collisionVec.Mag(),
                                                    0));
                } catch (std::runtime_error &e) {
                    std::cerr << "Unable to find correct un-centered distance (Mugast) : " << e.what() << std::endl;
                    localFragment.E_corrected[ii].push_back(tmpEn);
                    localFragment.BeamPosition_corrected[jj] = TVector3(0,0,0);
                }
            }
        }

        if ((reactionIt = reaction.find("M" + std::to_string(data->VAMOS_id_M) +
                                        "_Z" + std::to_string(data->VAMOS_id_Z) +
                                        "_" +
                                        localFragment.Particle[ii])) != reaction.end()){
            //Vamos and mugast localFragment found
            //Deformed Target
            localFragment.Ex[ii] = reactionIt->second->Set_Ek_Theta_Phi(localFragment.E[ii], 
                                                                        localFragment.EmissionDirection[ii].Theta(), 
                                                                        localFragment.EmissionDirection[ii].Phi());//WARNING this could be wring

            localFragment.E_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_Ek_cm();
            localFragment.Theta_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_P_cm().Vect().Theta();

            //Flat Target
            localFragment.ExFlatTarget[ii] = reactionIt->second->Set_Ek_Theta_Phi(localFragment.EFlatTarget[ii],
                                                                        localFragment.EmissionDirection[ii].Theta(),
                                                                        localFragment.EmissionDirection[ii].Phi());//WARNING this could be wring

            //Spherical Target
            localFragment.ExSphericalTarget[ii] = reactionIt->second->Set_Ek_Theta_Phi(localFragment.ESphericalTarget[ii],
                                                                        localFragment.EmissionDirection[ii].Theta(),
                                                                        localFragment.EmissionDirection[ii].Phi());//WARNING this could be wring

            localFragment.Ex_uncentered[ii].reserve(beamPositions.size());
            localFragment.Ex_corrected[ii].reserve(focusScale.size());
            if (ii == 0) {
                for (unsigned int jj = 0; jj < beamPositions.size(); ++jj) {
                    localFragment.Ex_uncentered[ii].push_back(reactionIt->second->Set_Ek_Theta(localFragment.E_uncentered[ii][jj],
                                                                                               localFragment.EmissionDirection_uncentered[jj].Theta()));
                }
                for (unsigned int jj = 0; jj < focusScale.size(); ++jj) {
                    localFragment.Ex_corrected[ii].push_back(reactionIt->second->Set_Ek_Theta(localFragment.E_corrected[ii][jj],
                                                                                              (localFragment.EmissionDirection_corrected[jj]).Angle(localFragment.BeamDirection_corrected[jj])));
                }
            }

        }
        else{
            if ((reactionIt = reaction.find(localFragment.Particle[ii])) != reaction.end()) {
                //Only mugast localFragment found
                localFragment.Ex[ii] = reactionIt->second->Set_Ek_Theta_Phi(localFragment.E[ii],
                                                                            localFragment.EmissionDirection[ii].Theta(),
                                                                            localFragment.EmissionDirection[ii].Phi());//WARNING this could be wring
                localFragment.E_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_Ek_cm();
                localFragment.Theta_CM[ii] = reactionIt->second->GetReactionFragment(3).Get_P_cm().Vect().Theta();

            }else
                localFragment.Ex[ii] = -2000;
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
                                                  200,                                //max_steps
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

double MugastIdentification::InitialBeamEnergy(double beamEnergyFromBrho, double tck){
    DEBUG("------------>Started: MugastIdentification::InitialBeamEnergy(), thickness : ", tck);
    beamEnergyFromBrho =
            energyLoss["beam"]["ice_back"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            tck * icePercentageSecond,
                                            0);

    beamEnergyFromBrho =
            energyLoss["beam"]["havar_back"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            havarThickness,
                                            0);

    beamEnergyFromBrho =
            energyLoss["beam"]["he3_front"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            averageBeamThickness,
                                            0);

    beamEnergyFromBrho =
            energyLoss["beam"]["havar_front"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            havarThickness,
                                            0);

    beamEnergyFromBrho =
            energyLoss["beam"]["ice_front"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            tck,
                                            0);
    DEBUG("------------>Finished: MugastIdentification::InitialBeamEnergy()", "");
    return beamEnergyFromBrho;
}

double MugastIdentification::middleTargetBeamEnergy(double beamEnergyFromBrho){
    beamEnergyFromBrho =
            energyLoss["beam"]["ice_back"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            currentIceThickness.first,
                                            0);

    beamEnergyFromBrho =
            energyLoss["beam"]["havar_back"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            havarThickness,
                                            0);

    beamEnergyFromBrho =
            energyLoss["beam"]["he3_front"]
                    ->EvaluateInitialEnergy(beamEnergyFromBrho,
                                            averageBeamThickness / 2.,
                                            0);

    return beamEnergyFromBrho;
}

TVector3 MugastIdentification::computeCollision(TVector3 emissionDirection, const TVector3& beamPosition) {
    double tmpThreashold{5E-6 * UNITS::mm};
    emissionDirection.SetMag(1);

    DEBUG("ComputeCollision","")
    DEBUG("beamPosition mag:", beamPosition.Mag())
    DEBUG("beamPosition theta:", beamPosition.Theta())
    DEBUG("beamPosition phi:", beamPosition.Phi())
    DEBUG("emissionDirection mag:", emissionDirection.Mag())
    DEBUG("emissionDirection theta:", emissionDirection.Theta())
    DEBUG("emissionDirection phi:", emissionDirection.Phi())

    //try{
    distanceMinimizer.reset(new Minimizer([&, this](double x) {
                                              emissionDirection.SetMag(x);
                                              auto vec = beamPosition + emissionDirection;
                                              if (emissionDirection.Z()>0)
                                                return pow(vec.Z() / UNITS::mm - this->gasThicknessCartesian->Evaluate(
                                                      sqrt(pow(vec.X() / UNITS::mm, 2) + pow(vec.Y() / UNITS::mm, 2))), 2);
                                              else
                                                return pow(vec.Z() / UNITS::mm + this->gasThicknessCartesian->Evaluate(
                                                      sqrt(pow(vec.X() / UNITS::mm, 2) + pow(vec.Y() / UNITS::mm, 2))), 2);
                                          },
                                          1.5 * UNITS::mm / abs(emissionDirection.CosTheta()),   //starting value
                                          3E-4 * UNITS::mm,                  //learning rate
                                          tmpThreashold,                     //threashold
                                          200,                               //max_steps
                                          0.98,                              //quenching
                                          1E-2 * UNITS::mm)                  //h
    );

    distanceMinimizer->Minimize();
    if (emissionDirection.Mag() > 15*UNITS::mm || isnan(emissionDirection.X())) {
        std::cerr   << "Found huge distance!!!!  " << std::endl
                    << "beam Mag : " << beamPosition.Mag()
                    << "at beam phi : " << beamPosition.Phi()<< ", theta : " << beamPosition.Theta() << std::endl;
        std::cerr   << "and emission vector phi : " << emissionDirection.Phi() << ", theta : " << emissionDirection.Theta()  << std::endl;

        throw std::runtime_error("Distance too high\n");
    }
    DEBUG("finished: ComputeCollision","")
    //}catch(...){
    //    return TVector3();
    //}
    return emissionDirection;
}

double MugastIdentification::computeAngleOfIncidence(const TVector3 &collisionVec, const TVector3 &beamPosition) {
    TVector3 vec = beamPosition + collisionVec;
    DEBUG("Angle of Incidence", "")
    DEBUG("collisionVec mag:", collisionVec.Mag())
    DEBUG("collisionVec theta:", collisionVec.Theta())
    DEBUG("collisionVec phi:", collisionVec.Phi())

    double x = sqrt(pow(vec.X() / UNITS::mm, 2) + pow(vec.Y() / UNITS::mm, 2));
    double derivative = gasThicknessDerivativeCartesian->Evaluate(x);//There should be no need to scale with UNITS::mm
    double angleFromDirection{0};
    double angleFromDerivative{0};
    if (collisionVec.Z() > 0){
        angleFromDirection = collisionVec.Angle(TVector3(0, 0, 1));
    }else{
        angleFromDirection = collisionVec.Angle(TVector3(0, 0, -1));
    }
    angleFromDerivative = atan(derivative);

    return angleFromDirection+angleFromDerivative;
}
