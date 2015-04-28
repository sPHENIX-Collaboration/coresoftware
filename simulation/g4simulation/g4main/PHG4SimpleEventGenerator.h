#ifndef PHG4SimpleEventGenerator_H__
#define PHG4SimpleEventGenerator_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <map>
#include <vector>

class TRandom3;
class PHG4InEvent;

class PHG4SimpleEventGenerator : public SubsysReco {

public:

  //! supported function distributions
  enum FUNCTION {Uniform,Gaus};

  PHG4SimpleEventGenerator(const std::string &name="EVTGENERATOR");
  virtual ~PHG4SimpleEventGenerator();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  //! random seed
  void set_seed(int seed);

  //! interface for adding particles by name
  void add_particles(std::string name, unsigned int count);

  //! interface for adding particle by pid
  void add_particles(int pid, unsigned int count);

  //! range of randomized eta values
  void set_eta_range(double eta_min, double eta_max);

  //! range of randomized phi values
  void set_phi_range(double phi_min, double phi_max);

  //! range of randomized pt values
  void set_pt_range(double pt_min, double mom_max);

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z);

  //! set the mean value of the vertex distribution
  void set_vertex_distribution_mean(double x, double y, double z);

  //! set the width of the vertex distribution function about the mean
  void set_vertex_distribution_width(double x, double y, double z);

  //! reuse the first existing vertex found
  void set_reuse_existing_vertex(bool b) {_reuse_existing_vertex = b;}

  //! set an offset vector from the existing vertex
  void set_existing_vertex_offset_vector(double x, double y, double z);
  
  //! set the distribution function of particles about the vertex
  void set_vertex_size_function(FUNCTION r);

  //! set the dimensions of the distribution of particles about the vertex
  void set_vertex_size_parameters(double mean, double width);

  //! raise the embed flag on the generated particles
  void set_embedflag(int embedflag) {_embedflag = embedflag;}

  //! print verbosity
  void set_verbosity(int verb) {Verbosity(verb);}

private:

  int get_pdgcode(std::string name);
  std::string get_pdgname(int pdgcode);
  double get_mass(int pid);

  int _seed;
  TRandom3 *_rand;

  // these need to be stored separately until run time when the names
  // can be translated using the GEANT4 lookup
  std::vector<std::pair<int, unsigned int> > _particle_codes; // <pdgcode, count>
  std::vector<std::pair<std::string, unsigned int> > _particle_names; // <names, count>
  bool _reuse_existing_vertex;
  FUNCTION _vertex_func_x;
  FUNCTION _vertex_func_y;
  FUNCTION _vertex_func_z;
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
  int _embedflag;

  PHG4InEvent* _ineve;
};

#endif
