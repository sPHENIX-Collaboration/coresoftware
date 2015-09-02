#ifndef PHG4ParticleGeneratorVectorMeson_H__
#define PHG4ParticleGeneratorVectorMeson_H__

#include "PHG4ParticleGeneratorBase.h"

class TRandom;
class TF1;

class PHG4ParticleGeneratorVectorMeson: public PHG4ParticleGeneratorBase
{
 public:
  PHG4ParticleGeneratorVectorMeson(const std::string &name="PGUN");
  virtual ~PHG4ParticleGeneratorVectorMeson() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  void set_eta_range(const double eta_min, const double eta_max);
  void set_rapidity_range(const double y_min, const double y_max);
  void set_mom_range(const double mom_min, const double mom_max);
  void set_pt_range(const double pt_min, const double pt_max);
  void set_vtx_zrange(const double zmin, const double zmax);
  void set_mass(const double mass);
  void set_width(const double width);
  void set_decay_types(const std::string &decay1, const std::string &decay2);
  void set_embedflag(int embedflag) {_embedflag = embedflag;}
  void set_histrand_init(const int initflag) {_histrand_init = initflag;}

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
  double width;
  double m1;
  double m2;
  int _embedflag;
  int _histrand_init;
  std::string decay1;
  std::string decay2;

  TF1 *fsin;
  TF1 *frap;
  TF1 *fpt;
  TRandom *trand;

};

#endif
