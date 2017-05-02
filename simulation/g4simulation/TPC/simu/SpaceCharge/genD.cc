#include <iostream>
#include <fstream>
#include <string>
#include "Langevin.h"

int main() {
  std::string kFileNameRoot="rho";
  float kInnerRadius = 30; //cm
  float kOutterRadius = 80; //cm
  float kHalfLength = 80; //cm
  int kNRadialSteps = 50;
  int kNAzimuthalSteps = 1;
  int kNLongitudinalSteps = 160;

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

  Langevin *langevin = new Langevin();
  langevin->SetDebugLevel(2);
  langevin->OutputFileName(kFileNameRoot);
  langevin->TPCDimensions( kInnerRadius, kOutterRadius, kHalfLength );
  langevin->TPCGridSize( kNRadialSteps, kNAzimuthalSteps, kNLongitudinalSteps );
  langevin->SetMirrorZ();
  langevin->Make();

  delete langevin;
  return 0;
}
