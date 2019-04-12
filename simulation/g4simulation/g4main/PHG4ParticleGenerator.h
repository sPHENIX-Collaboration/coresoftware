#ifndef PHG4ParticleGenerator_H__
#define PHG4ParticleGenerator_H__

#include "PHG4ParticleGeneratorBase.h"

class PHG4ParticleGenerator : public PHG4ParticleGeneratorBase
{
 public:
  PHG4ParticleGenerator(const std::string &name = "PGENERATOR");
  virtual ~PHG4ParticleGenerator() {}

  int process_event(PHCompositeNode *topNode);
  void set_z_range(const double z_min, const double z_max);
  void set_eta_range(const double eta_min, const double eta_max);
  void set_phi_range(const double phi_min, const double phi_max);
  void set_mom_range(const double mom_min, const double mom_max);
  void Print(const std::string &what = "ALL") const;

 protected:
  double z_min;
  double z_max;
  double eta_min;
  double eta_max;
  double phi_min;
  double phi_max;
  double mom_min;
  double mom_max;
};

#endif
