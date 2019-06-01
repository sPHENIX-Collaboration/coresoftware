// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SIMPLEEVENTGENERATOR_H
#define G4MAIN_PHG4SIMPLEEVENTGENERATOR_H

#include "PHG4ParticleGeneratorBase.h"

#include <map>
#include <vector>

class PHG4InEvent;
class PHCompositeNode;

class PHG4SimpleEventGenerator : public PHG4ParticleGeneratorBase {

public:

  //! supported function distributions
  enum FUNCTION {Uniform,Gaus};

  PHG4SimpleEventGenerator(const std::string &name="EVTGENERATOR");
  virtual ~PHG4SimpleEventGenerator(){}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  //! interface for adding particles by name
  void add_particles(const std::string &name, const unsigned int count);

  //! interface for adding particle by pid
  void add_particles(const int pid, const unsigned int count);

  //! set the starting time for the event
  void set_t0(const double t0);
  
  //! range of randomized eta values
  void set_eta_range(const double eta_min, const double eta_max);

  //! range of randomized phi values
  void set_phi_range(const double phi_min, const double phi_max);

  //! range of randomized pt values
  //! \param[in] pt_gaus_width   if non-zero, further apply a Gauss smearing to the pt_min - pt_max flat distribution
  void set_pt_range(const double pt_min, const double pt_max, const double pt_gaus_width = 0);

  //! range of randomized p values
  //! \param[in] p_gaus_width   if non-zero, further apply a Gauss smearing to the p_min - p_max flat distribution
  void set_p_range(const double p_min, const double p_max, const double p_gaus_width = 0);

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

private:

  double smearvtx(const double position, const double width, FUNCTION dist) const;
  // these need to be stored separately until run time when the names
  // can be translated using the GEANT4 lookup
  std::vector<std::pair<int, unsigned int> > _particle_codes; // <pdgcode, count>
  std::vector<std::pair<std::string, unsigned int> > _particle_names; // <names, count>  
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
  double _eta_min;
  double _eta_max;
  double _phi_min;
  double _phi_max;
  double _pt_min;
  double _pt_max;
  double _pt_gaus_width;
  double _p_min;
  double _p_max; 
  double _p_gaus_width;

  PHG4InEvent* _ineve;
};

#endif
