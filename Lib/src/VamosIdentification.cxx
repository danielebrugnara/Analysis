#include "VamosIdentification.h"

#ifdef VERBOSE_DEBUG
#  define DEBUG(x, y) std::cout << (x) << (y) << std::endl;
#else
#  define DEBUG(x, y) 
#endif

VamosIdentification::VamosIdentification() : cutsZ({18, 19, -1}),
                                             //cutsM({45, 46, 47, -1, -2}),
                                             //cuts_Q({13, 14, 15, 16, 17, 18, 19, -1, -2}),
                                             cutsM({45, 46, 47, -2}),
                                             //cuts_Q({13, 14, 15, 16, 17, 18, 19, -2}),
                                             data(nullptr)
{
}

bool VamosIdentification::initialize()
{
    DEBUG("------------>VamosIdentification::initialize()", "");

    //Focal plane aligments
    readFpTimeShifts();
    auto *tmpFile = new TFile("./Configs/Interpolations/InterpolationTimeFP.root");
    fpTimeInterpolation = new Interpolation(tmpFile);
    tmpFile->Close();

    //Masses in MeV [M][Z] as ints
    for (const auto &it_M : cutsM)
    {
        for (const auto &it_Z : cutsZ)
        {
            //mass[it_M][it_Z] = NPL::Nucleus(it_Z, it_M).Mass();
            if (it_Z > 0 && it_M > 0)
                mass[it_M][it_Z] = ReactionFragment(ReactionFragment::FragmentSettings(it_M, it_Z, it_Z, 0, 0)).Get_M();
            else
                mass[it_M][it_Z] = 0;
        }
    }

    auto *tmp = new std::unordered_map<std::string, TCutG *>();

    //It would be possible to parse the string to make automatic
    //what follows, however this gives possibility for better fine
    //tuning

    //Delta E vs Energy cuts

    DEBUG("Starting to look for VAMOS files", "");

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

    cut_type["dE2_E"] = *tmp;

    //Mass over Q cuts
    tmp = new std::unordered_map<std::string, TCutG *>();
    try
    {
        (*tmp)["MQ_Q_M46_Q18"] = cuts.at("MQ_Q_M46_Q18");
        (*tmp)["MQ_Q_M47_Q18"] = cuts.at("MQ_Q_M47_Q18");
        (*tmp)["MQ_Q_M45_Q16"] = cuts.at("MQ_Q_M45_Q16");
        (*tmp)["MQ_Q_M46_Q17"] = cuts.at("MQ_Q_M46_Q17");
        (*tmp)["MQ_Q_M47_Q17"] = cuts.at("MQ_Q_M47_Q17");
        (*tmp)["MQ_Q_M45_Q15"] = cuts.at("MQ_Q_M45_Q15");
        (*tmp)["MQ_Q_M45_Q17"] = cuts.at("MQ_Q_M45_Q17");
        (*tmp)["MQ_Q_M46_Q19"] = cuts.at("MQ_Q_M46_Q19");
        (*tmp)["MQ_Q_M46_Q16"] = cuts.at("MQ_Q_M46_Q16");
        (*tmp)["MQ_Q_M47_Q16"] = cuts.at("MQ_Q_M47_Q16");
        (*tmp)["MQ_Q_M45_Q14"] = cuts.at("MQ_Q_M45_Q14");
        (*tmp)["MQ_Q_M46_Q15"] = cuts.at("MQ_Q_M46_Q15");
        (*tmp)["MQ_Q_M47_Q15"] = cuts.at("MQ_Q_M47_Q15");
        (*tmp)["MQ_Q_M45_Q13"] = cuts.at("MQ_Q_M45_Q13");
        (*tmp)["MQ_Q_M46_Q14"] = cuts.at("MQ_Q_M46_Q14");
        (*tmp)["MQ_Q_M47_Q14"] = cuts.at("MQ_Q_M47_Q14");
        (*tmp)["MQ_Q_M45_Q12"] = cuts.at("MQ_Q_M45_Q12");
        (*tmp)["MQ_Q_M45_Q18"] = cuts.at("MQ_Q_M45_Q18");
        (*tmp)["MQ_Q_M47_Q19"] = cuts.at("MQ_Q_M47_Q19");
        //(*tmp)["MQ_Q_M-1_Q-1"] = cuts.at("MQ_Q_M-1_Q-1");
        (*tmp)["MQ_Q_M-2_Q-2"] = cuts.at("MQ_Q_M-2_Q-2");
    }
    catch (const std::out_of_range &err)
    {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the MQ_Q cuts\n");
    }

    cut_type["MQ_Q"] = *tmp;

    DEBUG("Found all needed cuts", "");

    return true;
}

VamosIdentification::~VamosIdentification()
{
    delete data;
    delete fpTimeInterpolation;
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

bool VamosIdentification::identify()
{
    fragment.BRho = **data->Brho;
    fragment.PTPosition.SetXYZ(**data->Pf,**data->Tf, 0);
    fragment.FocalPlanePosition.SetXYZ(**data->Xf,**data->Yf, 0);
    fragment.EmissionVersor.SetMagThetaPhi(1, **data->ThetaL, **data->PhiL);

    fragment.En = ((*data->IC)[0] > icThreashold) * ((*data->IC)[0] +
                                                     ((*data->IC)[1] > icThreashold) * (((*data->IC)[1] * (**data->Xf <= 45)) + (((*data->IC)[1] + 1.) * (**data->Xf > 45)) + //Correction for mis aligment IC1:Xf
                                                                                           ((*data->IC)[2] > icThreashold) * ((*data->IC)[2] +
                                                                                                                              ((*data->IC)[3] > icThreashold) * ((*data->IC)[3] +
                                                                                                                                                                 ((*data->IC)[4] > icThreashold) * ((*data->IC)[4] +
                                                                                                                                                                                                    ((*data->IC)[5] > icThreashold) * ((*data->IC)[5]))))));
    fragment.D_En = ((*data->IC)[0] > icThreashold) * ((*data->IC)[0] + ((*data->IC)[1] > icThreashold) * ((*data->IC)[1]));
    fragment.D_En2 = (*data->IC)[0] * ((*data->IC)[1] > icThreashold);

    //Computing the basic identifiaction
    fragment.T = getFpTime();
    fragment.Path = **data->Path + 5;
    fragment.V = fragment.Path / fragment.T;
    fragment.Beta = fragment.V / 29.9792;
    fragment.Gamma = 1. / sqrt(1.0 - fragment.Beta * fragment.Beta);
    fragment.M = (fragment.En) / 931.5016 / (fragment.Gamma - 1.);
    //mM2               = 18./20.8*(mE2)/931.5016/(mGamma2-1.);
    fragment.M_Q = **data->Brho / 3.105 / fragment.Beta / fragment.Gamma;
    fragment.Charge = fragment.M / fragment.M_Q;

    //dE - E identification
    for (const auto &z_search : cut_type.at("dE2_E"))
    {
        if (z_search.second->IsInside(fragment.En, fragment.D_En2))
        {
            //Z format dE2_E_Z18
            if (fragment.id_Z == 0)
                fragment.id_Z = std::stoi(
                    z_search.first.substr(z_search.first.find_last_of("_Z") + 1));
            else
                throw std::runtime_error("Overlapping Z gates\n");
        }
    }
    if (fragment.id_Z == 0)
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
    if (fragment.id_M == 0 || fragment.id_Q == 0)
        return false;

    //Lorentzvector computation
    fragment.p4.SetT(mass[fragment.id_M][fragment.id_Z]);
    TVector3 v4(0, 0, fragment.Beta);
    v4.SetMagThetaPhi(fragment.Beta, **data->ThetaL, **data->PhiL);
    fragment.p4.Boost(v4);

    return fragment.Identified = true;
}

double VamosIdentification::getShift()
{
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
    return shift + fpTimeInterpolation->Evaluate(**data->Xf);
}

double VamosIdentification::getEnFromBrho(){
    if(fragment.id_Z == 0 || fragment.id_M == 0 || fragment.id_Q == 0)
        return 0;
    return sqrt(pow(**data->Brho / 3.3356E-3 * fragment.id_Q, 2) + pow(mass[fragment.id_M][fragment.id_Z], 2)) - (mass[fragment.id_M][fragment.id_Q]);
}
