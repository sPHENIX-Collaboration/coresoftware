#include <iostream>
#include <fstream>
#include <string>
#include "QPileUpToy.h"

int main() {
  std::string kFileNameRoot="rho";
  float kInnerRadius = 20; //cm
  float kOutterRadius = 78; //cm
  float kHalfLength = 105.5; //cm
  int kNRadialSteps = 232;
  int kNAzimuthalSteps = 1;
  int kNLongitudinalSteps = 211;

  std::cout << "fetching geo.dat..." << std::endl;
  std::ifstream ifile("geo.dat");
  ifile >> kFileNameRoot;
  ifile >> kInnerRadius;  //cm
  ifile >> kOutterRadius; //cm
  ifile >> kHalfLength;   //cm
  ifile >> kNRadialSteps;
  ifile >> kNAzimuthalSteps;
  ifile >> kNLongitudinalSteps;
  ifile.close();
  std::cout << "Cage InnerRadius " << kInnerRadius << std::endl;
  std::cout << "Cage OuterRadius " << kOutterRadius << std::endl;
  std::cout << "Cage Longitudinal Length " << 2*kHalfLength << std::endl;
  std::cout << "==> ENLARGING granularity by 10" << std::endl;
  kNRadialSteps *= 10;
  kNAzimuthalSteps *= 10;
  kNLongitudinalSteps *= 10;
  std::cout << "Radial Steps " << kNRadialSteps << std::endl;
  std::cout << "Azimuthal Steps " << kNAzimuthalSteps << std::endl;
  std::cout << "Longitudinal Steps " << kNLongitudinalSteps << std::endl;
  ifile.close();

  float kGasFactor;
  float kMultiplicity;
  float kRate;
  float kIBF;
  std::cout << "fetching gas.dat..." << std::endl;
  ifile.open("gas.dat");
  ifile >> kGasFactor;
  ifile >> kMultiplicity;
  ifile >> kRate;
  ifile >> kIBF;
  ifile.close();

  QPileUp *initialdensity;
  initialdensity = new QPileUpToy(kGasFactor/*[Vs]*/, kMultiplicity, kRate/*[Hz]*/, kIBF); // PHENIX
  //initialdensity = new QPileUpToy(1.0/76628.0/*gasfactor [Vs]*/, 425.0 /*multiplicity*/, 5e+4/*rate [Hz]*/, 6); // PHENIX
  //initialdensity = new QPileUpToy(1.0/76628.0/*gasfactor [Vs]*/, 950.0 /*multiplicity*/, 5e+4/*rate [Hz]*/, 20, 1.5); // ALICE
  //initialdensity = new QPileUpToy(1.0/76628.0/*gasfactor [Vs]*/, 425.0 /*multiplicity*/, 15e+3/*rate [Hz]*/, 0); // STAR
  initialdensity->SetDebugLevel(2);
  initialdensity->OutputFileName(kFileNameRoot);
  initialdensity->TPCDimensions( kInnerRadius, kOutterRadius, kHalfLength );
  initialdensity->TPCGridSize( kNRadialSteps, kNAzimuthalSteps, kNLongitudinalSteps );
  initialdensity->Make();

  delete initialdensity;
  return 0;
}
