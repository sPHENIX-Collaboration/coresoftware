#include <iostream>
#include <fstream>
#include <string>
#include <TMath.h>
#include "FieldMaps.h"
#include "FieldMapsLaplace.h"

int main(int nvar, char ** cvar) {
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

  int binR = -1;
  TString par;
  if(nvar>1) {
    par = cvar[1];
    binR = par.Atoi();
    if(TMath::IsNaN(binR)) binR = -1;
  }
  std::cout << "routine for binR = " << binR << std::endl;

  FieldMaps *map;
  map = new FieldMapsLaplace();
  map->SetDebugLevel(1);
  map->OutputFileName(kFileNameRoot);
  map->TPCDimensions( kInnerRadius, kOutterRadius, kHalfLength );
  map->TPCGridSize( kNRadialSteps, kNAzimuthalSteps, kNLongitudinalSteps );
  map->MirrorZ();
  map->Make( binR );

  delete map;
  return 0;
}
