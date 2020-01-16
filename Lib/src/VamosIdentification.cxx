#include "VamosIdentification.h"

VamosIdentification::VamosIdentification() : cuts_Z({18, 19, -1}),
                                             //cuts_M({45, 46, 47, -1, -2}),
                                             //cuts_Q({13, 14, 15, 16, 17, 18, 19, -1, -2}), 
                                             cuts_M({45, 46, 47, -2}),
                                             cuts_Q({13, 14, 15, 16, 17, 18, 19, -2}), 
                                             data(nullptr),
                                             fragment(nullptr){}

bool VamosIdentification::Initialize() {
    ReadFPTimeShifts();

    //Masses in MeV [M][Z] as ints
    for (const auto &it_M : cuts_M) {
        for (const auto &it_Z : cuts_Z) {
            mass[it_M][it_Z] = it_M * AMU_TO_MEV;
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

    for (const auto &c : cuts) {
        std::cout << c.first << std::endl;
    }

#ifdef VERBOSE_DEBUG
    std::cout << "Starting to look for VAMOS files\n";
#endif

    try {
        (*tmp)["dE2_E_Z19"] = cuts.at("dE2_E_Z19");
        (*tmp)["dE2_E_Z18"] = cuts.at("dE2_E_Z18");
        (*tmp)["dE2_E_Z-1"] = cuts.at("dE2_E_Z-1");
    } catch (const std::out_of_range &err) {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the dE2_E cuts\n");
    }

    cut_type["dE2_E"] = *tmp;

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
        //(*tmp)["MQ_Q_M-1_Q-1"] = cuts.at("MQ_Q_M-1_Q-1");
        (*tmp)["MQ_Q_M-2_Q-2"] = cuts.at("MQ_Q_M-2_Q-2");
    } catch (const std::out_of_range &err) {
        std::cerr << err.what() << std::endl;
        throw std::runtime_error("Unable to find one of the MQ_Q cuts\n");
    }

    cut_type["MQ_Q"] = *tmp;

#ifdef VERBOSE_DEBUG
    std::cout << "Found all needed cuts\n";
#endif

    return true;
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
