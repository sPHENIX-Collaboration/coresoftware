// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SIMPLEEVENTGENERATOR_H
#define G4MAIN_PHG4SIMPLEEVENTGENERATOR_H

#include "PHG4ParticleGeneratorBase.h"

#include <cmath>
#include <map>
#include <string>   // for string
#include <utility>  // for pair
#include <vector>

class PHG4InEvent;
class PHCompositeNode;

class PHG4SimpleEventGenerator : public PHG4ParticleGeneratorBase
{
 public:
  //! supported function distributions
  enum FUNCTION
  {
    Uniform,
    Gaus
  };

  PHG4SimpleEventGenerator(const std::string &name = "EVTGENERATOR");
  ~PHG4SimpleEventGenerator() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  //! interface for adding particles by name
  void add_particles(const std::string &name, const unsigned int count);

  //! interface for adding particle by pid
  void add_particles(const int pid, const unsigned int count);

  //! range of randomized eta values (mutually exclusive with theta range)
  void set_eta_range(const double eta_min, const double eta_max);

  //! range of randomized theta values (mutually exclusive with eta range)
  void set_theta_range(const double theta_min, const double theta_max);

  //! range of randomized phi values
  void set_phi_range(const double phi_min, const double phi_max);

  //! power law value of distribution to sample from for pt values
  void set_power_law_n(const double n);

  //! range of randomized pt values (mutually exclusive with momentum range)
  //! \param[in] pt_gaus_width   if non-zero, further apply a Gauss smearing to the pt_min - pt_max flat distribution
  void set_pt_range(const double pt_min, const double pt_max, const double pt_gaus_width = 0);

  //! range of randomized p values (mutually exclusive with pt range)
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
  void set_vertex_size_function(FUNCTION r) { m_VertexSizeFunc_r = r; }

  //! set the dimensions of the distribution of particles about the vertex
  void set_vertex_size_parameters(const double mean, const double width);

 private:
  double smearvtx(const double position, const double width, FUNCTION dist) const;

  // these need to be stored separately until run time when the names
  // can be translated using the GEANT4 lookup
  std::vector<std::pair<int, unsigned int> > _particle_codes;          // <pdgcode, count>
  std::vector<std::pair<std::string, unsigned int> > _particle_names;  // <names, count>
                                                                       // so we can print out the function names without many if's
                                                                       // also used to check if function is implemented
  std::map<FUNCTION, std::string> m_FunctionNames = {{Uniform, "Uniform"}, {Gaus, "Gaus"}};

  PHG4InEvent *m_InEvent = nullptr;
  FUNCTION m_VertexFunc_x = Uniform;
  FUNCTION m_VertexFunc_y = Uniform;
  FUNCTION m_VertexFunc_z = Uniform;
  double m_Vertex_x = 0.;
  double m_Vertex_y = 0.;
  double m_Vertex_z = 0.;
  double m_VertexWidth_x = 0.;
  double m_VertexWidth_y = 0.;
  double m_VertexWidth_z = 0.;
  double m_VertexOffset_x = 0.;
  double m_VertexOffset_y = 0.;
  double m_VertexOffset_z = 0.;
  FUNCTION m_VertexSizeFunc_r = Uniform;
  double m_VertexSizeMean = 0.;
  double m_VertexSizeWidth = 0.;
  double m_EtaMin = -1.25;
  double m_EtaMax = 1.25;
  double m_ThetaMin = NAN;
  double m_ThetaMax = NAN;
  double m_PhiMin = -M_PI;
  double m_PhiMax = M_PI;
  double m_Pt_Min = 0.;
  double m_Pt_Max = 10.;
  double m_Pt_GausWidth = 0.;
  double m_P_Min = NAN;
  double m_P_Max = NAN;
  double m_P_GausWidth = NAN;
  double m_powerLawN = NAN;
};

#endif
