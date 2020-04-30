#ifndef __MINIMIZER_H__
#define __MINIMIZER_H__

#include <array>
#include <iostream>
#include <math.h>

class Minimizer
{
public:
  Minimizer(double (*) (const double &), double, double, double, int, double, double);
  ~Minimizer();

  //inline double PerformStep(double );
  double GetCoefficient() { return coefficient; };

private:
  double (*function_ptr) (const double &);
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
  double h;
  int n_steps;

public:
  double Minimize();
};
#endif
