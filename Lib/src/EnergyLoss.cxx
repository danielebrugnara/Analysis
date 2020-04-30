#include "EnergyLoss.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EnergyLoss::EnergyLoss()
{
  fInter = NULL;
  fSpline = NULL;
  fMax = -10000;
  fMin = 1e12;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EnergyLoss::~EnergyLoss()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EnergyLoss::EnergyLoss(std::string Path, std::string Source, int NumberOfSlice = 100, int LiseColumn, int NumberOfMass)
{
  fMax = -10000;
  fMin = 1e12;

  fNumberOfSlice = NumberOfSlice;
  fNumberOfMass = NumberOfMass;

#ifdef VERBOSE_DEBUG
  std::cout << "----->Initializing an EnergyLoss object " << std::endl;
#endif

  std::ifstream TableFile;

  // Opening dE/dX file
  TableFile.open(Path.c_str());
  if (TableFile.is_open())
  {
#ifdef VERBOSE_DEBUG
    std::cout << "Reading Energy Loss File : " << Path << std::endl;
#endif
  }
  else
  {
    std::cout << "ERROR: TABLE FILE " << Path << " NOT FOUND" << std::endl;
    exit(1);
  }

  if (Source == "G4Table" || Source == "G4table")
  {
    // Reading Data
    double energy, total;
    std::string dummy;
    //skipped first line
    getline(TableFile, dummy);
    while (TableFile >> energy >> total)
    {
      fEnergy.push_back(energy * UNITS::MeV);
      if (energy * UNITS::MeV > fMax)
      {
        fMax = energy * UNITS::MeV;
      }
      if (energy * UNITS::MeV < fMin)
      {
        fMin = energy * UNITS::MeV;
      }
      fdEdX_Total.push_back(total * (UNITS::MeV / UNITS::micrometer));
    }

    // Close File
    TableFile.close();
  }

  else if (Source == "SRIM")
  {
    // Reading Data
    double energy, nuclear, electronic;
    std::string unit, dummy;
    std::string line;
    while (std::getline(TableFile, line))
    {
      std::istringstream str(line);
      str >> energy >> unit >> electronic >> nuclear >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;
      if (energy == 0 || electronic == 0)
        continue;

      if (unit == "keV")
        energy = energy * UNITS::keV;
      if (unit == "MeV")
        energy = energy * UNITS::MeV;
      if (unit == "GeV")
        energy = energy * UNITS::GeV;
      fEnergy.push_back(energy);

      if (energy * UNITS::MeV > fMax)
        fMax = energy * UNITS::MeV;
      else if (energy * UNITS::MeV < fMin || fMin < 0)
        fMin = energy * UNITS::MeV;

      fdEdX_Nuclear.push_back(nuclear * UNITS::MeV / UNITS::mm);
      fdEdX_Electronic.push_back(electronic * UNITS::MeV / UNITS::mm);
      fdEdX_Total.push_back(nuclear * UNITS::MeV / UNITS::mm + electronic * UNITS::MeV / UNITS::mm);
    }

    // Close File
    TableFile.close();
  }

  else if (Source == "LISE")
  {
    // Reading Data
    double energy = 0, energyloss = 0;
    std::string dummy;
    // skipping comment first line
    getline(TableFile, dummy);

    while (TableFile >> energy)
    {
      for (int k = 0; k < 11; k++)
      {
        TableFile >> dummy;
        if (k + 1 == LiseColumn)
          energyloss = atof(dummy.c_str());
      }
      fEnergy.push_back(energy * UNITS::MeV);

      if (energy * UNITS::MeV > fMax)
        fMax = energy * UNITS::MeV;
      else if (energy * UNITS::MeV < fMin || fMin < 0)
        fMin = energy * UNITS::MeV;

      fdEdX_Total.push_back(energyloss * UNITS::MeV / UNITS::micrometer);
    }

    // Close File
    TableFile.close();
  }

  else
  {
    throw std::runtime_error("ERROR : EnergyLoss Wrong Source Type");
  }

  fInter = new TGraph(fEnergy.size(), &fEnergy[0], &fdEdX_Total[0]);

  fInter->Sort();

  fSpline = new TSpline3("Energy Loss Spline", fInter->GetX(), fInter->GetY(), fInter->GetN());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnergyLoss::Draw(std::string option) const
{
  fInter->Draw(option.c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::EvaluateNuclearLoss(double Energy) const
{
  if (fEnergy.size() == 0 || fdEdX_Nuclear.size() == 0)
  {
    std::cout << "No Nuclear table for this Energy loss";
    return -1000;
  }

  return Eval(Energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::EvaluateElectronicLoss(double Energy) const
{
  if (fEnergy.size() == 0 || fdEdX_Electronic.size() == 0)
  {
    std::cout << "No Electronic table for this Energy loss";
    return -1000;
  }

  return Eval(Energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::EvaluateTotalLoss(double Energy) const
{
  if (fEnergy.size() == 0 || fdEdX_Total.size() == 0)
  {
    std::cout << "No Total table for this Energy loss";
    return -1000;
  }

  return Eval(Energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EnergyLoss::Print() const
{
  std::cout << "Total Energy Loss : " << std::endl;
  int size = fdEdX_Total.size();
  for (int i = 0; i < size; i++)
    std::cout << fEnergy[i] / UNITS::MeV << " MeV "
              << fdEdX_Total[i] / UNITS::MeV * UNITS::micrometer << " MeV/um "
              << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::EnergyLossCalculation(double Energy,          // Energy of the detected particle
                                         double TargetThickness, // Target Thickness at 0 degree
                                         double Angle)           // Particle Angle
    const
{
  return (Energy - Slow(Energy, TargetThickness, Angle));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::Slow(double Energy,          // Energy of the detected particle
                        double TargetThickness, // Target Thickness at 0 degree
                        double Angle)           // Particle Angle
    const
{
  //   Lise file are given in MeV/u
  //   For SRIM and geant4 file fNumberOfMass = 1 whatever is the nucleus, file are given in MeV
  Energy = Energy / (double)fNumberOfMass;

  if (Angle > UNITS::CONSTANTS::halfpi)
    Angle = UNITS::CONSTANTS::pi - Angle;
  TargetThickness = TargetThickness / (cos(Angle));

  double SliceThickness = TargetThickness / (double)fNumberOfSlice;

  for (int i = 0; i < fNumberOfSlice; i++)
  {
    double de = Eval(Energy) * SliceThickness;
    Energy -= de / fNumberOfMass;

    if (Energy < 0)
    {
      Energy = 0;
      break;
    }
  }

  return (Energy * fNumberOfMass);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::EvaluateInitialEnergy(double Energy,          // Energy of the detected particle
                                         double TargetThickness, // Target Thickness at 0 degree
                                         double Angle)           // Particle Angle
    const
{
  if (TargetThickness == 0)
    return Energy;
  //   Lise file are given in MeV/u
  //   For SRIM and geant4 file fNumberOfMass = 1 whatever is the nucleus, file are given in MeV
  Energy = Energy / (double)fNumberOfMass;
  if (Angle > UNITS::CONSTANTS::halfpi)
    Angle = UNITS::CONSTANTS::pi - Angle;
  TargetThickness = TargetThickness / (cos(Angle));
  double SliceThickness = TargetThickness / (double)fNumberOfSlice;
  for (int i = 0; i < fNumberOfSlice; i++)
  {
    double de = Eval(Energy) * SliceThickness;
    Energy += de / fNumberOfMass;
  }
  return (Energy * fNumberOfMass);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::EvaluateMaterialThickness(double InitialEnergy, // Energy of the detected particle
                                             double FinalEnergy,
                                             double ThicknessLimit,
                                             double ThicknessStep) // Target Thickness at 0 degree
    const
{
  double Thickness = ThicknessStep;
  double Energy = InitialEnergy;
  while (Energy < FinalEnergy && Thickness > 0)
  {
    Energy = EvaluateInitialEnergy(Energy, ThicknessStep, 0);
    Thickness += ThicknessStep;
    if (Thickness > ThicknessLimit)
      Thickness = -1;
  }

  return Thickness;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//   Evaluate Total Energy of particle from Energy loss in a giver thickness
double EnergyLoss::EvaluateEnergyFromDeltaE(double DeltaE,           // Energy of the detected particle
                                            double TargetThickness,  // Target Thickness at 0 degree
                                            double Angle,            // Particle Angle
                                            double EnergyMin,        // Starting Energy
                                            double EnergyMax,        // Maximum Energy allowed
                                            double EnergyResolution, // Resolution at which function stop
                                            int MaxStep)             // Stop after MaxStep Whatever Precision is reached
    const
{

  double step_size = 10. * UNITS::MeV;
  double Energy = EnergyMax;
  double DE = 0;
  bool check_low = false;
  bool check_high = false;

  for (int i = 0; i < MaxStep; i++)
  {
    if (Energy > 0)
      DE = Energy - Slow(Energy, TargetThickness, Angle);
    else
      return 0;

    if (abs(DeltaE - DE) < EnergyResolution)
      return Energy;
    else if (DeltaE - DE > 0)
    {
      if (Energy - step_size > EnergyMin)
      {
        Energy = Energy - step_size;
        check_low = true;
      }
      else
      {
        step_size = step_size / 10.;
        Energy = Energy - step_size;
      }
    }

    else if (DeltaE - DE < 0)
    {
      if (Energy + step_size < EnergyMax)
      {
        Energy = Energy + step_size;
        check_high = true;
      }
      else
      {
        step_size = step_size / 10.;
        Energy = Energy + step_size;
      }
    }

    if (check_high && check_low)
    {
      step_size = step_size / 10.;
      check_high = false;
      check_low = false;
    }

    if (step_size < EnergyResolution)
      return Energy;
  }

  return Energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double EnergyLoss::Eval(double energy) const
{
  if (energy < 0)
  {
    throw std::runtime_error("\033[1;31mWARNING: negative energy given to EnergyLoss\033[0m\n");
    return 0;
  }
  else if (energy != energy)
  {
    throw std::runtime_error("WARNING: nan energy given to EnergyLoss\n");
    return 0;
  }

  return fInter->Eval(energy, fSpline);
}
