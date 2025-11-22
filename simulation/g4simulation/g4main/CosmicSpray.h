// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_COSMICSPRAY_H
#define G4MAIN_COSMICSPRAY_H

#include "EcoMug.h"

#include <fun4all/SubsysReco.h>

#include <limits>
#include <string>  // for string

class PHCompositeNode;

// class CosmicSpray : public PHG4ParticleGeneratorBase
class CosmicSpray : public SubsysReco
{
 public:
  bool InDetector(double x, double y, double z) const;
  CosmicSpray(const std::string &name = "COSMICS", const double R = 650);
  ~CosmicSpray() override = default;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  void set_gen_min_momentum(const double p) { gen.SetMinimumMomentum(p); }
  void set_gen_max_momentum(const double p) { gen.SetMaximumMomentum(p); }

 private:
  EcoMug gen;

  double _gun_e = std::numeric_limits<double>::quiet_NaN();
  //  double _x_min = std::numeric_limits<double>::quiet_NaN();
  double _x_max = std::numeric_limits<double>::quiet_NaN();
  //  double _z_min = std::numeric_limits<double>::quiet_NaN();
  double _z_max = std::numeric_limits<double>::quiet_NaN();
  double _y_fix = std::numeric_limits<double>::quiet_NaN();

  double _R = std::numeric_limits<double>::quiet_NaN();
};
#endif
