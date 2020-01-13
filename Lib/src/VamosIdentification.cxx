#include "VamosIdentification.h"

VamosIdentification::VamosIdentification() : cuts_Z({18, 19, -1}),
                                             cuts_M({45, 46, 47, -1, -2}),
                                             cuts_Q({13, 14, 15, 16, 17, 18, 19, -1, -2}) {
    ReadFPTimeShifts();

    //Masses in MeV [M][Z] as ints
    for (const auto &it_M : cuts_M) {
        for (const auto &it_Z : cuts_Z) {
            mass[it_M][it_Z] = it_M *AMU_TO_MEV;
        }
    }

    //More precise data
    mass[46][18] = 45.968082712 * AMU_TO_MEV;  //in MeV
    mass[47][18] = 46.972934865 * AMU_TO_MEV;  //in MeV
    mass[47][19] = 46.961661614 * AMU_TO_MEV;  //in MeV
    mass[46][19] = 45.961981586 * AMU_TO_MEV;  //In MeV

    std::unordered_map<std::string, TCutG *> *tmp = new std::unordered_map<std::string, TCutG *>();

    //It would be possible to parse the string to make automatic
    //what follows, however this gives possibility for better fine
    //tuning

    //Delta E vs Energy cuts
    try {
        (*tmp)["dE2_E_Z19"] = cuts.at("dE2_E_Z19");
        (*tmp)["dE2_E_Z18"] = cuts.at("dE2_E_Z18");
        (*tmp)["dE2_E_Z-1"] = cuts.at("dE2_E_Z-1");
    } catch (const std::out_of_range &err) {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the dE2_E cuts\n");
    }

    cut_type["DE2_E"] = *tmp;

    //Mass over Q cuts
    tmp = new std::unordered_map<std::string, TCutG *>();
    try {
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
        (*tmp)["MQ_Q_M-1_Q-1"] = cuts.at("MQ_Q_M-1_Q-1");
        (*tmp)["MQ_Q_M-2_Q-2"] = cuts.at("MQ_Q_M-2_Q-2");
    } catch (const std::out_of_range &err) {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the MQ_Q cuts\n");
    }

    cut_type["MQ_Q"] = *tmp;
}

VamosIdentification::~VamosIdentification() {}

void VamosIdentification::ReadFPTimeShifts() {
    std::ifstream cali_file("./Configs/FP_Time.cal");
    std::string line;
    std::string min;
    std::string max;
    std::string shift;
    if (!cali_file) throw std::runtime_error("Error Opening cali file\n");
    while (std::getline(cali_file, line)) {
        std::stringstream str;
        str << line;
        str >> min >> max >> shift;
        TimeShifts.push_back(std::pair<double, double>(stod(max), stod(shift)));
    }
    double tmp = -1E6;
    for (const auto &it : TimeShifts) {
        if (it.first < tmp) throw std::runtime_error("Time shifts not incremental\n");
        tmp = it.first;
    }
}

inline double VamosIdentification::GetShift() {
    double shift{0};
    for (const auto &it : TimeShifts) {
        if (**data->Xf < it.first)
            shift = it.first;
        else
            break;
    }
    return shift;
}

inline double VamosIdentification::GetFPTime() {
    //    return 540.5*(data->AGAVA_VAMOSTS<104753375647998)+537.9*(data->AGAVA_VAMOSTS>=104753375647998) -2.*data->T_FPMW_CATS2_C + GetShift();
    return 0;
}

inline bool VamosIdentification::Identify() {
    fragment->En = ((*data->IC)[0] > IC_threashold) * ((*data->IC)[0] + ((*data->IC)[1] > IC_threashold) * ((*data->IC)[1] + ((*data->IC)[2] > IC_threashold) * ((*data->IC)[2] + ((*data->IC)[3] > IC_threashold) * ((*data->IC)[3] + ((*data->IC)[4] > IC_threashold) * ((*data->IC)[4] + ((*data->IC)[5] > IC_threashold) * ((*data->IC)[5]))))));
    fragment->D_En = ((*data->IC)[0] > IC_threashold) * ((*data->IC)[0] + ((*data->IC)[1] > IC_threashold) * ((*data->IC)[1]));
    fragment->D_En2 = (*data->IC)[0] * ((*data->IC)[1] > IC_threashold);

    //fragment->T = 540.5 * (data->AGAVA_VAMOSTS < 104753375647998) + 537.9 * (data->AGAVA_VAMOSTS >= 104753375647998) - 2. * (data->T_FPMW_CATS2_C) + 2.7 * (MW_N[0] == 16) + 2.7 * (MW_N[0] == 15) + 2.9 * (MW_N[0] == 14) + 2.9 * (MW_N[0] == 13) + 2.4 * (MW_N[0] == 12) + 1.3 * (MW_N[0] == 11) + 1.5 * (MW_N[0] == 10) + 1.6 * (MW_N[0] == 9) - 0.6 * (MW_N[0] == 8) + 2.5 * (MW_N[0] == 7) + 2. * (MW_N[0] == 6) + 1.6 * (MW_N[0] == 5) + 1.1 * (MW_N[0] == 4) - 0.6 * (MW_N[0] == 3) - 1.2 * (MW_N[0] == 2) - 4.0 * (MW_N[0] == 1);
    fragment->T = GetFPTime();
    fragment->Path = **data->Path + 5;

    //Computing the basic identifiaction
    fragment->V = fragment->Path / fragment->T;
    fragment->Beta = fragment->V / 29.9792;
    fragment->Gamma = 1. / sqrt(1.0 - fragment->Beta * fragment->Beta);
    fragment->M = (fragment->En) / 931.5016 / (fragment->Gamma - 1.);
    //mM2 = 18./20.8*(mE2)/931.5016/(mGamma2-1.);
    fragment->M_Q = **data->Brho / 3.105 / fragment->Beta / fragment->Gamma;
    fragment->Charge = fragment->M / fragment->M_Q;

    //dE - E identification
    for (const auto &z_search : cut_type.at("dE2_E")) {
        if (z_search.second->IsInside(fragment->En, fragment->D_En2)) {
            //Z format dE2_E_Z18
            if (fragment->id_Z == 0)
                fragment->id_Z = stoi(
                    z_search.first.substr(z_search.first.find_last_of("_Z") + 1));
            else
                throw std::runtime_error("Overlapping Z gates\n");
        }
    }
    if (fragment->id_Z == 0) return false;

    //MQ - Q identification
    for (const auto &mq_search : cut_type.at("MQ_Q")) {
        if (mq_search.second->IsInside(fragment->M_Q, fragment->Charge)) {
            if (fragment->id_M == 0 && fragment->id_Q == 0) {
                fragment->id_M = stoi(mq_search.first.substr(mq_search.first.find_last_of("M") + 1, 2));
                fragment->id_Q = stoi(mq_search.first.substr(mq_search.first.find_last_of("_") + 2));
            } else
                throw std::runtime_error("Overlapping M_Q gates");
        }
    }
    if (fragment->id_M == 0 || fragment->id_Q == 0) return false;

    //Lorentzvector computation
    fragment->p4.SetT(mass[fragment->id_M][fragment->id_Z]);
    TVector3 v4(0, 0, fragment->Beta);
    v4.SetMagThetaPhi(fragment->Beta, **data->ThetaL, **data->PhiL);
    fragment->p4.Boost(v4);

    return fragment->Identified = true;
}

inline void VamosIdentification::SetData(Data const *data) {
    delete this->data;
    delete this->fragment;
    fragment = new Fragment();
    this->data = data;
}
