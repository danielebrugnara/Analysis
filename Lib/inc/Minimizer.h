#ifndef __MINIMIZER_H__
#define __MINIMIZER_H__

#include <array>
#include <iostream>
#include <math.h>

class Minimizer
{
public:
  Minimizer(double, double, double, double, double, int, double);
  ~Minimizer();

  //inline double PerformStep(double );
  double GetCoefficient() { return coefficient; };

private:
  std::array<double, 2> y;
  std::array<double, 2> x;
  std::array<long double, 2> derivative;
  std::array<long double, 2> step;
  std::array<long double, 2> rate;
  double startingRate;
  double coefficient;
  double threshold;
  double max_steps;
  double quenching;
  int n_steps;

public:
  inline double PerformStep(double new_value)
  {
    new_value = pow(new_value, 2);
    if (new_value == 0)
      return x[0];
    if (abs(new_value) < threshold)
    { //Already ok
#ifdef VERBOSE_DEBUG
      std::cout << "Value within threashold : " << new_value << std::endl;
#endif
      return x[1];
    }
    y[0] = y[1];
    y[1] = new_value;
    rate[0] = rate[1];
    step[0] = step[1];
    derivative[0] = derivative[1];
    if (n_steps > 1)
    {
      if (y[1] - y[0] == 0)
        derivative[1] = derivative[0] * 0.5;
      else
        derivative[1] = (y[1] - y[0]) / (step[0]);
      rate[1] = rate[0];
      if (abs(derivative[1]) < 1.E-3)
        rate[1] = rate[1] * 2;
      step[1] = (-1.) * quenching * derivative[1] * rate[1];
      //std::cout << step[1] <<std::endl;
      //if (derivative[1]<1.E-7 && derivative[0]>1.E-3)
      //    step[1] = 2* step[0];
      //std::cout << step[1] <<std::endl;
      //if (abs(step[1])> 5 )
      //    step[1] = ((step[1] > 0) - (step[1] < 0)) * step[0] * 2;
      //std::cout << step[1] <<std::endl;
      quenching *= quenching;
    };
    n_steps++;
    x[0] = x[1];
    x[1] = x[0] + step[1];
#ifdef VERBOSE_DEBUG
    std::cout << n_steps
              << "------> Derivative : " << derivative[1]
              << "\t\t uno : " << y[0] << "\t\t due : " << y[1]
              //<<"\t\t difference : " << _X[1] - _X[0]
              //<<"\t\t old step : " << _Step[0]
              << "\t\t new step : " << step[1]
              << "\t\t position : " << x[1] + step[1]
              << std::endl;
    if (rate[0] != rate[1])
      std::cout << "Changed rate from :" << rate[0] << " To : " << rate[1] << std::endl;
#endif
    return x[1];
  }
};

#endif
