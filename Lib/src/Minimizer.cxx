#include "Minimizer.h"


Minimizer::Minimizer(double start_value, 
                        double start_coeff = 0.,  
                        double start_step = 0.1, 
                        double learning_rate = 1.5, 
                        double thresh = 0.1, 
                        int max_steps = 100, 
                        double quenching = 1):
                            y({0, start_value}),
                            x({0, start_coeff}),
                            derivative({0, 0}),
                            step({0, start_step}),
                            rate({learning_rate, learning_rate}),
                            startingRate(learning_rate),
                            threshold(thresh),
                            max_steps(max_steps),
                            quenching(quenching),
                            n_steps(0){
}

Minimizer::~Minimizer(){
}