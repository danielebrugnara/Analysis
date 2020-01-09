#include "VamosIdentification.h"

VamosIdentification::VamosIdentification(){
    ReadFPTimeShifts();
}

VamosIdentification::~VamosIdentification(){}

void VamosIdentification::ReadFPTimeShifts(){
	std::ifstream cali_file("./Configs/FP_Time.cal");
	std::string line;
	std::string min;
	std::string max;
	std::string shift;
	if (!cali_file) throw std::runtime_error("Error Opening cali file\n");
	while (std::getline(cali_file, line)){
		std::stringstream str;
		str << line;
		str >> min >> max >> shift;
        TimeShifts.push_back(std::pair<double, double>(stod(max), stod(shift)));
	}
    double tmp = -1E6;
    for(const auto & it: TimeShifts){
        if (it.first<tmp) throw std::runtime_error("Time shifts not incremental\n"); 
        tmp = it.first;
    }
}

inline double VamosIdentification::GetShift(){
    double shift {0};
    for(const auto & it: TimeShifts){
        if (data->Xf < it.first) shift = it.first;
        else break;
    }
    return shift;
}


inline double VamosIdentification::GetFPTime(){
    return 540.5*(data->AGAVA_VAMOSTS<104753375647998)+537.9*(data->AGAVA_VAMOSTS>=104753375647998) -2.*data->T_FPMW_CATS2_C + GetShift();
}


inline void VamosIdentification::Identify(){

    fragment->En = (data->IC[0]>IC_threashold)*( data->IC[0] + (data->IC[1]>IC_threashold)*(data->IC[1] +(data->IC[2]>IC_threashold)*( data->IC[2] +(data->IC[3]>IC_threashold)*( data->IC[3] +(data->IC[4]>IC_threashold)*( data->IC[4] +(data->IC[5]>IC_threashold)*( data->IC[5]))))));
    fragment->D_En = (data->IC[0]>IC_threashold)*( data->IC[0] + (data->IC[1]>IC_threashold)*( data->IC[1]));
    fragment->D_En2 = data->IC[0] *(data->IC[1]> IC_threashold);

    //fragment->T = 540.5 * (data->AGAVA_VAMOSTS < 104753375647998) + 537.9 * (data->AGAVA_VAMOSTS >= 104753375647998) - 2. * (data->T_FPMW_CATS2_C) + 2.7 * (MW_N[0] == 16) + 2.7 * (MW_N[0] == 15) + 2.9 * (MW_N[0] == 14) + 2.9 * (MW_N[0] == 13) + 2.4 * (MW_N[0] == 12) + 1.3 * (MW_N[0] == 11) + 1.5 * (MW_N[0] == 10) + 1.6 * (MW_N[0] == 9) - 0.6 * (MW_N[0] == 8) + 2.5 * (MW_N[0] == 7) + 2. * (MW_N[0] == 6) + 1.6 * (MW_N[0] == 5) + 1.1 * (MW_N[0] == 4) - 0.6 * (MW_N[0] == 3) - 1.2 * (MW_N[0] == 2) - 4.0 * (MW_N[0] == 1);
    fragment->T = GetFPTime();
    fragment->Path = data->Path + 5;

    //Computing the basic identifiaction
    fragment->V = fragment->Path / fragment->T;
    fragment->Beta = fragment->V / 29.9792;
    fragment->Gamma = 1. / sqrt(1.0 - fragment->Beta * fragment->Beta);
    fragment->M = (fragment->En) / 931.5016 / (fragment->Gamma - 1.);
    //mM2 = 18./20.8*(mE2)/931.5016/(mGamma2-1.);
    fragment->M_Q = data->Brho / 3.105 / fragment->Beta / fragment->Gamma;
    fragment->Charge = fragment->M / fragment->M_Q;
}
