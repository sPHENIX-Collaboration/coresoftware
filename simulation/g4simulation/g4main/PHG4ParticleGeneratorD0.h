// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEGENERATORD0_H
#define G4MAIN_PHG4PARTICLEGENERATORD0_H

#include "PHG4ParticleGeneratorBase.h"

#include <string>                       // for string

class PHCompositeNode;
class TF1;

class PHG4ParticleGeneratorD0 : public PHG4ParticleGeneratorBase
{
 public:
  PHG4ParticleGeneratorD0(const std::string &name = "D0GEN");
  virtual ~PHG4ParticleGeneratorD0() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  void set_eta_range(const double eta_min, const double eta_max);
  void set_rapidity_range(const double y_min, const double y_max);
  void set_mom_range(const double mom_min, const double mom_max);
  void set_pt_range(const double pt_min, const double pt_max);
  void set_vtx_zrange(const double zmin, const double zmax);
  void set_mass(const double mass);

 protected:
  double vtx_zmin;
  double vtx_zmax;
  double y_min;
  double y_max;
  double eta_min;
  double eta_max;
  double mom_min;
  double mom_max;
  double pt_min;
  double pt_max;
  double mass;
  double m1;
  double m2;

  TF1 *fsin;
  TF1 *frap;
  TF1 *fpt;
};

#endif
