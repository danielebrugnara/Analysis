#include "Minimizer.h"

Minimizer::Minimizer(double (*function_ptr)(const double &),
                     double starting_value,
                     double learning_rate = 0.7,
                     double threshold = 0.1,
                     int max_steps = 100,
                     double quenching = 1,
                     double h=1E-3)  :      
                            function_ptr(function_ptr),
                            y({0, 0}),
                            x({0, starting_value}),
                            derivative({0, 0}),
                            step({0, 0}),
                            rate({learning_rate, learning_rate}),
                            threshold(threshold),
                            max_steps(max_steps),
                            quenching(quenching),
                            h(h),
                            n_steps(0)
{
}

Minimizer::~Minimizer()
{
}


double Minimizer::Minimize()
{
	double fx ;
	do {
        derivative[0]=derivative[1];
		derivative[1]=(function_ptr(x[0]+h)-function_ptr(x[0]-h))/(2*h);
        step[0]=step[1];
		step[1]=-derivative[1]*rate[1];
        fx = function_ptr(x[0]);
		do {
			x[1] = x[0]+step[1];
			//if(x[1]<a) {
			//	x[1]=a;
			//} else if(x1>b) {
			//	x[1]=b;
			//}
			if(function_ptr(x[1]) >fx) {
				step[1]/=2;
				continue;
			}
            if(derivative[1] < 1E-3){
                step[1]*=2;
                continue;
            }
		} while(0);
        rate[0]=rate[1];
		rate[1]=rate[1]*quenching;
        if(++n_steps>max_steps)
            throw std::runtime_error("Reached max steps in minimizer\n");
	} while(fabs(x[1]-x[0])>threshold);
    return x[1];
}
