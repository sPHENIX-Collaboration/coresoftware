// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PARTICLEGENERATORVECTORMESON_H
#define G4MAIN_PHG4PARTICLEGENERATORVECTORMESON_H

#include "PHG4ParticleGeneratorBase.h"

#include <cmath>
#include <map>
#include <string>  // for string

class PHCompositeNode;
class PHG4InEvent;
class TRandom;
class TF1;

class PHG4ParticleGeneratorVectorMeson : public PHG4ParticleGeneratorBase
{
 public:
  //! supported function distributions
  enum FUNCTION
  {
    Uniform,
    Gaus
  };

  explicit PHG4ParticleGeneratorVectorMeson(const std::string &name = "VMESON");
  ~PHG4ParticleGeneratorVectorMeson() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  //! interface for adding particles by name
  void add_decay_particles(const std::string &name1, const std::string &name2, const unsigned int decay_id);
  void add_decay_particles(const std::string &name, const unsigned int decay_id);

  void set_decay_vertex_offset(double dx, double dy, double dz, const unsigned int decay_id);
  void set_eta_range(const double eta_min, const double eta_max);
  void set_rapidity_range(const double y_min, const double y_max);
  void set_mom_range(const double mom_min, const double mom_max);
  void set_pt_range(const double pt_min, const double pt_max);
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

  void set_read_vtx_from_hepmc(bool read_vtx) { read_vtx_from_hepmc = read_vtx; }

  void set_mass(const double mass_in) { mass = mass_in; }
  void set_width(const double width_in) { m_Width = width_in; }
  void set_decay_types(const std::string &decay1, const std::string &decay2);
  void set_histrand_init(const int initflag) { _histrand_init = initflag; }
  void set_upsilon_1s();
  void set_upsilon_2s();
  void set_upsilon_3s();

 private:
  double smearvtx(const double position, const double width, FUNCTION dist) const;
  std::map<unsigned int, int> decay1_codes;          // <pdgcode, count>
  std::map<unsigned int, std::string> decay1_names;  // <names, count>
  std::map<unsigned int, int> decay2_codes;          // <pdgcode, count>
  std::map<unsigned int, std::string> decay2_names;  // <names, count>
  std::map<unsigned int, double> decay_vtx_offset_x;
  std::map<unsigned int, double> decay_vtx_offset_y;
  std::map<unsigned int, double> decay_vtx_offset_z;

  FUNCTION _vertex_func_x = Uniform;
  FUNCTION _vertex_func_y = Uniform;
  FUNCTION _vertex_func_z = Uniform;
  double _vertex_x = 0.;         // primary vertex (or mean) x component, cf. vtx_x = track-by-track vertex x component
  double _vertex_y = 0.;         // primary vertex (or mean) y component, cf. vtx_y = track-by-track vertex y component
  double _vertex_z = 0.;         // primary vertex (or mean)z component, cf. vtx_z = track-by-track vertex z component
  double _vertex_width_x = 0.;   // sigma x if not use existing vtx
  double _vertex_width_y = 0.;   // sigma y if not use existing vtx
  double _vertex_width_z = 0.;   // sigma z if not use existing vtx
  double _vertex_offset_x = 0.;  // track-by-track decay vertex offset if use existing vtx
  double _vertex_offset_y = 0.;  // track-by-track decay vertex offset if use existing vtx
  double _vertex_offset_z = 0.;  // track-by-track decay vertex offset if use existing vtx
  FUNCTION _vertex_size_func_r = Uniform;
  double _vertex_size_mean = 0.;
  double _vertex_size_width = 0.;
  bool read_vtx_from_hepmc = true;

  double y_min = 0.;
  double y_max = 0.;
  double eta_min = -1.;
  double eta_max = 1.;
  double mom_min = 0.;
  double mom_max = 10.;
  double pt_min = 4.;
  double pt_max = 4.;
  double mass = NAN;
  double m_Width = NAN;
  double m1 = NAN;
  double m2 = NAN;
  int _histrand_init = 0;
  std::string decay1 = "e+";
  std::string decay2 = "e-";

  TF1 *fsin = nullptr;
  TF1 *frap = nullptr;
  TF1 *fpt = nullptr;
  TRandom *trand = nullptr;

  PHG4InEvent *ineve = nullptr;
};

#endif
