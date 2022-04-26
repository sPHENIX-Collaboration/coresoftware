// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_COSMICSPRAY_H
#define G4MAIN_COSMICSPRAY_H

#include <g4main/PHG4ParticleGeneratorBase.h>

#include <cmath>
#include <map>
#include <string>   // for string
#include <utility>  // for pair
#include <vector>
#include "TF3.h"
#include <g4main/PHG4InEvent.h>
#include "EcoMug.h"
class PHG4InEvent;
class PHCompositeNode;

class CosmicSpray : public PHG4ParticleGeneratorBase
{
public:
  bool InDetector(double x, double y, double z);
  CosmicSpray(const std::string &name, const double R, const int &debug);
  ~CosmicSpray() override {}
  int process_event(PHCompositeNode *topNode) override;
  
 private:

 EcoMug gen;
 
  static double _gun_e;
  static double _x_min;
  static double _x_max;
  static double _z_min;
  static double _z_max;
  static double _y_fix;
  
  int _debug;
  double _R;

  PHG4InEvent *_InEvent = nullptr;
 
};
#endif
