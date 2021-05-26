// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEGENERATOR_H
#define G4MAIN_PHG4PARTICLEGENERATOR_H

#include "PHG4ParticleGeneratorBase.h"

#include <cmath>
#include <string>  // for string

class PHCompositeNode;

class PHG4ParticleGenerator : public PHG4ParticleGeneratorBase
{
 public:
  PHG4ParticleGenerator(const std::string &name = "PGENERATOR");
  ~PHG4ParticleGenerator() override {}

  int process_event(PHCompositeNode *topNode) override;
  void set_z_range(const double z_min, const double z_max);
  void set_eta_range(const double eta_min, const double eta_max);
  void set_phi_range(const double phi_min, const double phi_max);
  void set_mom_range(const double mom_min, const double mom_max);
  void Print(const std::string &what = "ALL") const override;

 protected:
  double m_ZMin = -10.;
  double m_ZMax = 10.;
  double m_EtaMin = -1.;
  double m_EtaMax = 1.;
  double m_PhiMin = -M_PI;
  double m_PhiMax = M_PI;
  double m_MomMin = 0.;
  double m_MomMax = 10.;
};

#endif
