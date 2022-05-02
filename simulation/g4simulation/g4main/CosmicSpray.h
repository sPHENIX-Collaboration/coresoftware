// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_COSMICSPRAY_H
#define G4MAIN_COSMICSPRAY_H

#include "EcoMug.h"

#include <fun4all/SubsysReco.h>

#include <cmath>
#include <string>  // for string

class PHCompositeNode;

//class CosmicSpray : public PHG4ParticleGeneratorBase
class CosmicSpray : public SubsysReco
{
 public:
  bool InDetector(double x, double y, double z);
  CosmicSpray(const std::string &name = "COSMICS", const double R = 650);
  ~CosmicSpray() override {}
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

 private:
  EcoMug gen;

  double _gun_e = NAN;
  double _x_min = NAN;
  double _x_max = NAN;
  double _z_min = NAN;
  double _z_max = NAN;
  double _y_fix = NAN;

  double _R = NAN;
};
#endif
