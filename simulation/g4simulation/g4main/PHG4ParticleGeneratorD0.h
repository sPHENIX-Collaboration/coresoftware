// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEGENERATORD0_H
#define G4MAIN_PHG4PARTICLEGENERATORD0_H

#include "PHG4ParticleGeneratorBase.h"

#include <string>  // for string

class PHCompositeNode;
class TF1;

class PHG4ParticleGeneratorD0 : public PHG4ParticleGeneratorBase
{
 public:
  PHG4ParticleGeneratorD0(const std::string &name = "D0GEN");
  ~PHG4ParticleGeneratorD0() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_eta_range(const double eta_min, const double eta_max);
  void set_rapidity_range(const double y_min, const double y_max);
  void set_mom_range(const double mom_min, const double mom_max);
  void set_pt_range(const double pt_min, const double pt_max);
  void set_vtx_zrange(const double zmin, const double zmax);
  void set_mass(const double mass);

 protected:
  double vtx_zmin = -10.;
  double vtx_zmax = 10.;
  double y_min = 0.;
  double y_max = 0.;
  double eta_min = -1.;
  double eta_max = 1.;
  double mom_min = 0.;
  double mom_max = 10.;
  double pt_min = 4.;
  double pt_max = 4.;
  double mass = 1.86486;
  double m1 = 0.493677;
  double m2 = 0.13957018;

  TF1 *fsin = nullptr;
  TF1 *frap = nullptr;
  TF1 *fpt = nullptr;
};

#endif
