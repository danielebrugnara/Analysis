#include "VamosIdentification.h"

#ifdef VERBOSE_DEBUG
#  define DEBUG(x, y) std::cout << (x) << (y) << std::endl;
#else
#  define DEBUG(x, y) 
#endif

VamosIdentification::VamosIdentification() : cutsZ({18, 19, -1}),
                                             //cutsM({45, 46, 47, -1, -2}),
                                             //cuts_Q({13, 14, 15, 16, 17, 18, 19, -1, -2}),
                                             cutsM({44, 45, 46, 47, 48, -2}),
                                             //cuts_Q({13, 14, 15, 16, 17, 18, 19, -2}),
                                             data(nullptr)
{
}

bool VamosIdentification::initialize()
{
    DEBUG("------------>VamosIdentification::initialize()", "");

    //Focal plane aligments
    readFpTimeShifts();
    auto *tmpFPFile = new TFile("./Configs/Interpolations/InterpolationPosFP.root");
    fpPosInterpolation = new Interpolation(tmpFPFile);
    tmpFPFile->Close();

    auto *tmpQFile = new TFile("./Configs/Interpolations/M46_Time.root");
    massInterpolation = new Interpolation(tmpQFile, true);
    tmpQFile->Close();

    //Masses in MeV [M][Z] as ints
    for (const auto &it_M : cutsM)
    {
        for (const auto &it_Z : cutsZ)
        {
            //mass[it_M][it_Z] = NPL::Nucleus(it_Z, it_M).Mass();
            if (it_Z > 0 && it_M > 0)
                mass[it_M][it_Z] = new ReactionFragment(ReactionFragment::FragmentSettings(it_M, it_Z, it_Z, 0, 0));
            else
                mass[it_M][it_Z] = new ReactionFragment(ReactionFragment::FragmentSettings(1, 1, 1, 0, 0));
        }
    }


    //It would be possible to parse the string to make automatic
    //what follows, however this gives possibility for better fine
    //tuning

    //Delta E vs Energy cuts

    DEBUG("Starting to look for VAMOS files", "");

    //De-E 2
    auto *tmp = new std::unordered_map<std::string, TCutG *>();
    try
    {
        (*tmp)["dE2_E_Z19"] = cuts.at("dE2_E_Z19");
        (*tmp)["dE2_E_Z18"] = cuts.at("dE2_E_Z18");
        (*tmp)["dE2_E_Z-1"] = cuts.at("dE2_E_Z-1");
    }
    catch (const std::out_of_range &err)
    {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the dE2_E cuts\n");
    }
    cut_type["dE2_E"] = *tmp;

    //De-E 1
    delete tmp;
    tmp = new std::unordered_map<std::string, TCutG *>();
    try
    {
        (*tmp)["dE_E_Z19"] = cuts.at("dE_E_Z19");
        (*tmp)["dE_E_Z18"] = cuts.at("dE_E_Z18");
        (*tmp)["dE_E_Z-1"] = cuts.at("dE_E_Z-1");
    }
    catch (const std::out_of_range &err)
    {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the dE_E cuts\n");
    }
    cut_type["dE_E"] = *tmp;


    //Mass over Q cuts
    delete tmp;
    tmp = new std::unordered_map<std::string, TCutG *>();
    try{
        for (const auto& it: cuts){
            if (it.first.compare(0, 5, "MQ_Q_") != 0) continue;
            (*tmp)[it.first] = it.second;
        }
    }catch (const std::out_of_range &err){
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the MQ_Q cuts\n");
    }

    cut_type["MQ_Q"] = *tmp;
    delete tmp;

    DEBUG("Found all needed cuts", "");

    //Ionization chamber calibration
    icCalibration[0][0] = 0;
    icCalibration[1][0] = -2;
    icCalibration[2][0] = -6.4;
    icCalibration[3][0] = -5;
    icCalibration[4][0] = -4.5;
    icCalibration[5][0] = 0;
    icCalibration[6][0] = 0;

    icCalibration[0][1] = 1;
    icCalibration[1][1] = 1.02;
    icCalibration[2][1] = 1.072;
    icCalibration[3][1] = 1.04;
    icCalibration[4][1] = 0.8;
    icCalibration[5][1] = 1;
    icCalibration[6][1] = 1;

    return true;
}

VamosIdentification::~VamosIdentification()
{
    delete data;
    delete fpPosInterpolation;
    delete massInterpolation;

    for(auto& it1: mass){
        for(auto& it2: it1.second)
            delete it2.second;
    }
}

void VamosIdentification::readFpTimeShifts()
{
    std::ifstream caliFile("./Configs/Calibrations/FP_Time.cal");
    std::string line;
    std::string min;
    std::string max;
    std::string shift;
    if (!caliFile)
        throw std::runtime_error("Error Opening cali file\n");
    while (std::getline(caliFile, line))
    {
        std::stringstream str;
        str << line;
        str >> min >> max >> shift;
        timeShifts.emplace_back(std::stod(max), std::stod(shift));
    }
    double tmp = -1E6;
    DEBUG("Dumping FP shift calibration", "");
    for (const auto &it : timeShifts)
    {
        if (it.first < tmp)
            throw std::runtime_error("Time shifts not incremental\n");
        tmp = it.first;
        DEBUG("it.first = ", it.first );
        DEBUG("it.second = ", it.second );
    }
}

bool VamosIdentification::identify(){
    fragment.BRho = **data->Brho;
    fragment.PTPosition.SetXYZ(**data->Pf,**data->Tf, 0);
    fragment.FocalPlanePosition.SetXYZ(**data->Xf,**data->Yf, 0);
    fragment.EmissionVersor.SetMagThetaPhi(1, **data->ThetaL, **data->PhiL);

    for(int i{0}; i<7; ++i)
        fragment.IC[i] = ((*data->IC)[i])*icCalibration[i][1]+icCalibration[i][0];

    //fragment.En = ((*data->IC)[0] > icThreashold) * ((*data->IC)[0] +
    //                                                 ((*data->IC)[1] > icThreashold) * (((*data->IC)[1] * (**data->Xf <= 45)) + (((*data->IC)[1] + 1.) * (**data->Xf > 45)) + //Correction for mis aligment IC1:Xf
    //                                                                                       ((*data->IC)[2] > icThreashold) * ((*data->IC)[2] +
    //                                                                                                                          ((*data->IC)[3] > icThreashold) * ((*data->IC)[3] +
    //                                                                                                                                                             ((*data->IC)[4] > icThreashold) * ((*data->IC)[4] +
    //                                                                                                                                                                                                ((*data->IC)[5] > icThreashold) * ((*data->IC)[5]))))));

    fragment.En = 0;
    for(int i{0}; i<6 && fragment.IC[i]>icThreashold; ++i)
        fragment.En += fragment.IC[i];

    for(int i{0}; i<2 && fragment.IC[i]>icThreashold; ++i)
        fragment.D_En += fragment.IC[i];

    for(int i{0}; i<1 && fragment.IC[i]>icThreashold; ++i)
        fragment.D_En2 += fragment.IC[i];
    //fragment.D_En = ((*data->IC)[0] > icThreashold) * ((*data->IC)[0] + ((*data->IC)[1] > icThreashold) * ((*data->IC)[1]));
    //fragment.D_En2 = (*data->IC)[0] * ((*data->IC)[1] > icThreashold);

    //Computing the basic identifiaction
    fragment.T = getFpTime();
    fragment.Path = **data->Path/cos(**data->Pf/1000.) + (760.-752.81)/(cos(**data->Pf/1000.)*cos(**data->Tf/1000.))+fpPosInterpolation->Evaluate(**data->Xf);
    fragment.V = fragment.Path / fragment.T;
    fragment.Beta = fragment.V / 29.9792;
    fragment.Gamma = 1. / sqrt(1.0 - fragment.Beta * fragment.Beta);
    fragment.En_NM = 1.2902850 *fragment.IC[0] + 11.441050; //From fit of simulation
    fragment.M_Q = **data->Brho / 3.105 / fragment.Beta / fragment.Gamma;
    double massCorrectionFactor = 46-massInterpolation->Evaluate(**data->TW);
    //double correctionFactor0 = 14.3656; //From comparing energy with brhoenergy
    //double correctionFactor1 = 0.902245;
    double correctionFactor0 = 12.6699; //From comparing energy with brhoenergy
    double correctionFactor1 = 0.924928;
    //fragment.M =  (correctionFactor1*(fragment.En + fragment.En_NM)+correctionFactor0) / 931.5016 / (fragment.Gamma - 1.);
    fragment.M = massCorrectionFactor + (correctionFactor1*(fragment.En + fragment.En_NM)+correctionFactor0) / 931.494 / (fragment.Gamma - 1.);

    //Correction on M_Q vs M graph
    correctionFactor0 = 42.3969;
    correctionFactor1 = 1.34352;
    fragment.M = fragment.M - (fragment.M_Q*correctionFactor1+correctionFactor0) +46;

    //mM2               = 18./20.8*(mE2)/931.5016/(mGamma2-1.);
    fragment.Charge = fragment.M / fragment.M_Q;

    //dE2 - E identification
    for (const auto &z_search : cut_type.at("dE2_E"))
    {
        if (z_search.second->IsInside(fragment.En, fragment.D_En2))
        {
            //Z format dE2_E_Z18
            if (fragment.id2_Z == 0)
                fragment.id2_Z = std::stoi(
                    z_search.first.substr(z_search.first.find_last_of("_Z") + 1));
            else
                throw std::runtime_error("Overlapping Z gates\n");
        }
    }

    //dE - E identification
    for (const auto &z_search : cut_type.at("dE_E"))
    {
        if (z_search.second->IsInside(fragment.En, fragment.D_En))
        {
            //Z format dE_E_Z18
            if (fragment.id_Z == 0)
                fragment.id_Z = std::stoi(
                        z_search.first.substr(z_search.first.find_last_of("_Z") + 1));
            else
                throw std::runtime_error("Overlapping Z gates\n");
        }
    }

    if (fragment.id2_Z == 0 && fragment.id_Z == 0)
        return false;

    //MQ - Q identification
    for (const auto &mq_search : cut_type.at("MQ_Q"))
    {
        if (mq_search.second->IsInside(fragment.M_Q, fragment.Charge))
        {
            if (fragment.id_M == 0 && fragment.id_Q == 0)
            {
                fragment.id_M = std::stoi(mq_search.first.substr(mq_search.first.find_last_of("M") + 1, 2));
                fragment.id_Q = std::stoi(mq_search.first.substr(mq_search.first.find_last_of("_") + 2));
            }
            else
                throw std::runtime_error("Overlapping M_Q gates :  " + mq_search.first +
                                         "  and  M" + std::to_string(fragment.id_M) +
                                         "  Q" + std::to_string(fragment.id_Q) + "\n");
        }
    }
    fragment.En_BRho = mass[46][18]->GetEkFromBrho_Q(fragment.BRho, 17);
    if (fragment.id_M == 0 || fragment.id_Q == 0)
        return false;

    //Lorentzvector computation
    fragment.p4.SetT(mass[fragment.id_M][fragment.id2_Z !=0 ? fragment.id2_Z : fragment.id_Z]->Get_M());
    TVector3 v4(0, 0, fragment.Beta);
    v4.SetMagThetaPhi(fragment.Beta, **data->ThetaL, **data->PhiL);
    fragment.p4.Boost(v4);

    //En from brho computation
    //TODO: correct this!!!
    //fragment.En_BRho = mass[fragment.id_M][fragment.id2_Z !=0 ? fragment.id2_Z : fragment.id_Z]->GetEkFromBrho_Q(fragment.BRho, fragment.id_Q);

    return fragment.Identified = true;
}

double VamosIdentification::getShift(){
    double shift{timeShifts.at(0).second};
    if (**data->Xf == -1500 || **data->Xf > timeShifts.back().first)
        return 0;
    for (const auto &it : timeShifts)
    {
        if (**data->Xf > it.first)
            shift = it.second;
        else
        {
            shift = it.second;
            break;
        }
    }
    //return shift + fpPosInterpolation->Evaluate(**data->Xf);
    return shift;
}

double VamosIdentification::getFpTime(){
    return  540.5 * (**data->AGAVA_VAMOSTS < 104753375647998)
            + 537.9 * (**data->AGAVA_VAMOSTS >= 104753375647998) - 2. * **data->T_FPMW_CATS2_C
            +getShift();
}

//double VamosIdentification::getEnFromBrho(){
//    if(fragment.id_Z == 0 || fragment.id_M == 0 || fragment.id_Q == 0)
//        return 0;
//    return sqrt(pow(**data->Brho / 3.3356E-3 * fragment.id_Q, 2) + pow(mass[fragment.id_M][fragment.id_Z]->Get_M(), 2)) - (mass[fragment.id_M][fragment.id_Q]->Get_M());
//}
