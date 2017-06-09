#ifndef PHG4ParticleGeneratorVectorMeson_H__
#define PHG4ParticleGeneratorVectorMeson_H__

#include "PHG4ParticleGeneratorBase.h"

class TRandom;
class TF1;

class PHG4ParticleGeneratorVectorMeson: public PHG4ParticleGeneratorBase
{
 public:
  //! supported function distributions
      enum FUNCTION {Uniform,Gaus};

  PHG4ParticleGeneratorVectorMeson(const std::string &name="PGUN");
  virtual ~PHG4ParticleGeneratorVectorMeson() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  void set_eta_range(const double eta_min, const double eta_max);
  void set_rapidity_range(const double y_min, const double y_max);
  void set_mom_range(const double mom_min, const double mom_max);
  void set_pt_range(const double pt_min, const double pt_max);
  void set_vtx_zrange(const double zmin, const double zmax);
  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z);

  //! set the mean value of the vertex distribution
  void set_vertex_distribution_mean(const double x, const double y, const double z);

  //! set the width of the vertex distribution function about the mean
  void set_vertex_distribution_width(const double x, const double y, const double z);

  //! set an offset vector from the existing vertex
  void set_existing_vertex_offset_vector(const double x, const double y, const double z);

  //! set the distribution function of particles about the vertex
  void set_vertex_size_function(FUNCTION r);

  //! set the dimensions of the distribution of particles about the vertex
  void set_vertex_size_parameters(const double mean, const double width);

  void set_mass(const double mass);
  void set_width(const double width);
  void set_decay_types(const std::string &decay1, const std::string &decay2);
  void set_histrand_init(const int initflag) {_histrand_init = initflag;}

 private:

  double smearvtx(const double position, const double width, FUNCTION dist) const;

 protected:
  double vtx_zmin;
  double vtx_zmax;
  FUNCTION _vertex_func_x;
  FUNCTION _vertex_func_y;
  FUNCTION _vertex_func_z;
  double _t0;
  double _vertex_x;
  double _vertex_y;
  double _vertex_z;
  double _vertex_width_x;
  double _vertex_width_y;
  double _vertex_width_z;
  double _vertex_offset_x;
  double _vertex_offset_y;
  double _vertex_offset_z;
  FUNCTION _vertex_size_func_r;
  double _vertex_size_mean;
  double _vertex_size_width;

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
  int _histrand_init;
  std::string decay1;
  std::string decay2;

  TF1 *fsin;
  TF1 *frap;
  TF1 *fpt;
  TRandom *trand;

};

#endif
