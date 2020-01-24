#include "Calibration.h"

Calibration::Calibration(std::string file_name, 
                            unsigned int n_parameters,
                            unsigned int channels): 
                                n_parameters(n_parameters),
                                n_channels(channels){
    if (channels<1 || n_parameters<2) throw std::runtime_error("Calibration parameters incorrect");

    coefficients.resize(n_channels);
    for(unsigned int ii = 0; ii<n_channels; ++ii ){
        coefficients.at(ii).resize(n_parameters);
        for(unsigned int jj = 0; jj<n_parameters; ++jj ){
            coefficients.at(ii).at(jj) = 0;
        }
    }        

    std::ifstream calibration_file(file_name);
    if (!calibration_file.is_open()) calibration_is_present = false;
    std::string line;
    int channel;
    double coeff;
    while (std::getline(calibration_file, line)){
        std::istringstream str (line);
        str >> channel;
        if (channel >= n_channels) continue;
        int ii=0;
        while(str >> coeff){
            coefficients[channel][ii++] = coeff;
        }
    }
}

Calibration::~Calibration(){
}