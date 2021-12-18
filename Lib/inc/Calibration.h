#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <math.h>

class Calibration
{
public:
    Calibration(std::string, unsigned int, unsigned int, bool=true);
    ~Calibration();

    inline double Evaluate(const double &value, const int &channel){
        if (!channelMode) throw std::runtime_error("Something is wrong in calibration, this should not have been called\n");
        double result = value;
        for (long unsigned int ii = 0; ii < coefficients[channel].size(); ++ii){
            result = result + coefficients[channel][ii] * pow(value, ii);
        }
        return result;
    }
    inline double Evaluate(const double &value, const double &parameter){
        if (channelMode) throw std::runtime_error("Something is wrong in calibration, this should not have been called\n");
        unsigned int idx{0};
        for(unsigned int i{1}; i<intervals.size() && parameter>intervals[i-1]; ++i)
            idx = i;
        double result = value;
        for (long unsigned int ii = 0; ii < coefficients[idx].size(); ++ii){
            result = result + coefficients[idx][ii] * pow(value, ii);
        }
        return result;
    }

private:
    bool calibrationIsPresent;
    bool channelMode;
    const int nParameters;
    const unsigned int nChannels;
    std::vector<std::vector<double>> coefficients;
    std::vector<double> intervals;
};