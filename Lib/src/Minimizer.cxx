#include "Minimizer.h"


Minimizer::Minimizer(double start_value, bool verbose = false, double start_coeff = 0.,  
                    double start_step = 0.1, double learning_rate = 1.5, 
                     double thresh = 0.001, int max_steps = 100, 
                     double quenching = 1):
                        Y({0, start_value}),
                        X({0, start_coeff}),
                        Derivative({0, 0}),
                        Step({0, start_step}),
                        Rate({learning_rate, learning_rate}),
                        StartingRate(learning_rate),
                        Threshold(thresh),
                        MaxSteps(max_steps),
                        Quenching(quenching),
                        NSteps(0),
                        Verbose(verbose)
                    {
}