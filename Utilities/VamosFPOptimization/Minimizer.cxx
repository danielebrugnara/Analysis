#include "Minimizer.h"

//Minimizer::Minimizer(double (*function_ptr)(const double &),
Minimizer::Minimizer(std::function<double(const double &)> function_ptr,
                     double starting_value,
                     double learning_rate = 0.7,
                     double threshold = 0.1,
                     int max_steps = 100,
                     double quenching = 1,
                     double h=1E-3)  :      
                            function_ptr(function_ptr),
                            x({starting_value, starting_value}),
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

Minimizer::Minimizer(double starting_value,
                     double learning_rate = 0.7,
                     double threshold = 0.1,
                     int max_steps = 100,
                     double quenching = 1,
                     double h=1E-3)  :      
                            function_ptr(function_ptr),
                            x({starting_value, starting_value}),
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

void Minimizer::SetThreashold(double threashold){
    this->threshold=threashold;
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
        x[0] = x[1];
	    x[1] = x[0]+step[1];
#ifdef VERBOSE_DEBUG
        	std::cout   << "Step number : " << n_steps << "--->"
                	    << " x value : "<< x[0]
                	    << " y value : "<< fx
                	    << " derivative value : "<< derivative[1]
                	    << " step value : " << step[1]
                	    << std::endl;
#endif
		if(function_ptr(x[1]) >fx && n_steps>2) {
			step[1]/=2;
            break;
		}
        rate[0]=rate[1];
		rate[1]=rate[1]*quenching;
        ++n_steps;
        if(n_steps>max_steps)
            	throw std::runtime_error("Reached max steps in minimizer\n");
	} while(fabs(x[1]-x[0])>threshold);
    return x[1];
}


double Minimizer::Step(std::pair<double, double> y_value){
    derivative[0]=derivative[1];
	derivative[1]=(y_value.second-y_value.first)/h;
    step[0]=step[1];
	step[1]=-derivative[1]*rate[1];
    x[0] = x[1];
	x[1] = x[0]+step[1];

    rate[0]=rate[1];
	rate[1]=rate[1]*quenching;

    ++n_steps;
    return x[1];
} 

