#pragma once

//Eloss calculations from NPTOOL: https://gitlab.in2p3.fr/np/nptool

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <stdlib.h>

//   ROOT
#include "TGraph.h"
#include "TSpline.h"
#include "TAxis.h"
#include "TFile.h"

// Use CLHEP System of unit and Physical Constant
#include "Units.h"

class EnergyLoss
{

public: //   Constructor
  EnergyLoss();
  EnergyLoss(std::string Path,      //   Path of dE/dX table file
             std::string Source,    // Type of file : Geant4,Lise,SRIM
             int NumberOfSlice,     //   Low number = Faster, High Number = more accurate / typical: 100 to 1000
             int LiseColumns = 0,   //   Indicate which model to read in a lise File, set to 0 (Default value) for a SRIM / Geant4 file
             int NumberOfMass = 1); //   Number of mass A of the nucleus (used only for Lise file)
  ~EnergyLoss();

private:                                //   dE/dX, slice parameter
  int fNumberOfSlice;                   // Number Of Slice used to evaluate energy loss
  int fNumberOfMass;                    // Number of mass A of the nucleus (used only for Lise file)
  double fMax;                          // Max energy value of the table
  double fMin;                          // Min energy value of the tahble
  std::vector<double> fEnergy;          // Energy
  std::vector<double> fdEdX_Nuclear;    // Nuclear Stopping Power
  std::vector<double> fdEdX_Electronic; // Electronic Stopping Power
  std::vector<double> fdEdX_Total;      // Total Stopping Power
  //Interpolator*     fInter            ; // Interpolator Used to evaluate Energy loss at given energy
  TGraph *fInter;    // Graph use to perform interpolation
  TGraph *fInterNuclear;    // Graph use to perform interpolation
  TGraph *fInterElectronic;    // Graph use to perform interpolation
  TSpline3 *fSpline; // Spline 3rd order used to perfom the interpolation

public:                           //   General Function on dE/dX table
  double Eval(double ener) const; // return the evaluated energy
  double EvaluateNuclearLoss(double ener) const;
  double EvaluateElectronicLoss(double ener) const;
  double EvaluateTotalLoss(double ener) const;

public: //   Function to Slow down particle or reconstruct their initial energy
  //   Calculate Energy of a particle after crossing material
  double Slow(double Energy,          // Energy of the detected particle
              double TargetThickness, // Target Thickness at 0 degree
              double Angle)           // Particle Angle
      const;

  //   Calculate Energy Loss of a particle inside a material
  double EnergyLossCalculation(double Energy,          // Energy of the detected particle in MeV
                               double TargetThickness, // Target Thickness at 0 degree in mm
                               double Angle)           // Particle Angle in rad
      const;

  //   Evaluate Initial Energy of particle before crossing material knowing Angle, final Energy
  //   and Target Thickness.
  double EvaluateInitialEnergy(double energy,          // Energy of the detected particle
                               double TargetThickness, // Target Thickness at 0 degree
                               double Angle)           // Particle Angle
      const;

  //   Evaluate Total Energy of particle from Energy loss in a giver thickness
  double EvaluateEnergyFromDeltaE(double DeltaE,           // Energy of the detected particle
                                  double TargetThickness,  // Target Thickness at 0 degree
                                  double Angle,            // Particle Angle
                                  double EnergyMin,        // Starting Energy
                                  double EnergyMax,        // Maximum Energy allowed
                                  double EnergyResolution, // Resolution at which function stop
                                  int MaxStep = 1000000)   // Stop after MaxStep Whatever Precision is reached
      const;

  // Evaluate the thickness the particle has been through using the energy loss and initial energy
  // usefull for thickness measurement using particle sources
  double EvaluateMaterialThickness(double InitialEnergy, // Energy of the detected particle
                                   double FinalEnergy,
                                   double ThicknessLimit,
                                   double ThicknessStep = 0.1 * UNITS::micrometer) // Target Thickness at 0 degree
      const;

public: //   Usefull function
  //   Display parameter
  void Print() const;
  //   Draw (CERN ROOT)
  void Draw(const std::string& option = "") const;
  void DrawNuclear(const std::string& option = "") const;
  void DrawElectronic(const std::string& option = "") const;
};