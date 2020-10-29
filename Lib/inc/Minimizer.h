#pragma once

#include <array>
#include <functional>
#include <iostream>
#include <math.h>

class Minimizer
{
public:
  Minimizer(std::function<double(const double &)>, double, double, double, int, double, double);
  ~Minimizer();

  void SetThreashold(double);
private:
  //double (*function_ptr) (const double &);
  std::function<double(const double &)> function_ptr;
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