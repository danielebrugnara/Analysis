#ifndef __CALIBRATION_H__
#define __CALIBRATION_H__

#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include <math.h>

class Calibration {
   public:
    Calibration(std::string, unsigned int, unsigned int);
    ~Calibration();

    inline double Evaluate(const double & value, const int & channel){
        double result = 0;
        for (int ii=0; ii<coefficients[channel].size(); ++ii){
            result = result + coefficients[channel][ii] * pow(value, ii);
        }
        return result;
    }

   private:
    bool calibration_is_present;
    const int n_parameters;
    const int n_channels;
    std::vector<std::vector<double>> coefficients; 
};

#endif