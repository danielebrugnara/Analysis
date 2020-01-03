#include "Minimizer.h"
#include <array>

Minimizer::Minimizer(double start_value, bool verbose = false, double start_coeff = 0.,  
                    double start_step = 0.1, double learning_rate = 1.5, 
                     double thresh = 0.001, int max_steps = 100, 
                     double quenching = 1):
                        _Y({0, start_value}),
                        _X({0, start_coeff}),
                        _Derivative({0, 0}),
                        _Step({0, start_step}),
                        _Rate({learning_rate, learning_rate}),
                        _StartingRate(learning_rate),
                        _Threshold(thresh),
                        _MaxSteps(max_steps),
                        _Quenching(quenching),
                        _NSteps(0),
                        _Verbose(verbose)
                    {
}

Minimizer::~Minimizer(){
}

inline double Minimizer::Step(double new_value)
{
  new_value = abs(new_value);
  if (new_value == 0)
    return _X[0];
  if (abs(new_value) < _Threshold)
    return _X[1];
  _Y[0] = _Y[1];
  _Y[1] = new_value;
  _Rate[0] = _Rate[1];
  _Step[0] = _Step[1];
  _Derivative[0] = _Derivative[1];
  if (_NSteps > 1)
  {
    if (_Y[1] - _Y[0] == 0)
      _Derivative[1] = _Derivative[0] * 0.5;
    else
      _Derivative[1] = (_Y[1] - _Y[0]) / (_Step[0]);
    _Rate[1] = _Rate[0];
    if (abs(_Derivative[1])<1.E-3) _Rate[1] = _Rate[1]*2;
    _Step[1] = (-1.) * _Quenching * _Derivative[1] * _Rate[1];
    if (_Derivative[1]<1.E-7 && _Derivative[0]>1.E-3) _Step[1] = 2*_Step[0];
    if (abs(_Step[1])> 5 ) _Step[1] = _Step[0] * 2;
    _Quenching *= _Quenching;
    if (_Verbose)
    {
      std::cout << _NSteps
                << "------> Derivative : " << _Derivative[1]
                << "\t\t uno : " << _Y[0] << "\t\t due : " << _Y[1]
                //<<"\t\t difference : " << _X[1] - _X[0]
                //<<"\t\t old step : " << _Step[0]
                << "\t\t new step : " << _Step[1]
                << "\t\t position : " << _X[1] + _Step[1]
                << std::endl;
      if (_Rate[0]!=_Rate[1]) std::cout <<"Changed rate from :" <<_Rate[0] << " To : " << _Rate[1] <<std::endl;
    }
  };
  _NSteps++;
  _X[0] = _X[1];
  _X[1] = _X[0] + _Step[1];
  return _X[1];
}
