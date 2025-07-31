// $Id: $

/*!
 * \file PHHepMCGenHelper.h
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHHEPMC_PHHEPMCGENHELPER_H
#define PHHEPMC_PHHEPMCGENHELPER_H

#include <CLHEP/Vector/ThreeVector.h>

#include <gsl/gsl_rng.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

class PHCompositeNode;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;

namespace HepMC
{
  class GenEvent;
}

/*!
 * \brief PHHepMCGenHelper provides service of DST upload of HepMC subevent, vertex assignment and random generator
 */
class PHHepMCGenHelper
{
 public:
  PHHepMCGenHelper();
  virtual ~PHHepMCGenHelper();

  //! supported function distributions
  enum VTXFUNC
  {
    //! uniform distribution with half width set via set_vertex_distribution_width()
    Uniform,
    //! normal distribution with sigma width set via set_vertex_distribution_width()
    Gaus
  };

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(VTXFUNC x, VTXFUNC y, VTXFUNC z, VTXFUNC t);

  //! set the mean value of the vertex distribution, use PHENIX units of cm, ns
  void set_vertex_distribution_mean(const double x, const double y, const double z, const double t);

  //! set the width of the vertex distribution function about the mean, use PHENIX units of cm, ns
  void set_vertex_distribution_width(const double x, const double y, const double z, const double t);

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return _embedding_id; }
  //
  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { _embedding_id = id; }
  //
  //! reuse vertex from another PHHepMCGenEvent with embedding_id = src_embedding_id Additional smearing and shift possible with set_vertex_distribution_*()
  void set_reuse_vertex(int src_embedding_id)
  {
    _reuse_vertex = true;
    _reuse_vertex_embedding_id = src_embedding_id;
  }

  //
  //! init interface nodes
  virtual int create_node_tree(PHCompositeNode *topNode);

  //! choice of reference version of the PHHepMCGenEvent
  static const PHHepMCGenEvent *get_PHHepMCGenEvent_template();

  //! send HepMC::GenEvent to DST tree. This function takes ownership of evt
  PHHepMCGenEvent *insert_event(HepMC::GenEvent *evt);

  const PHHepMCGenEventMap *get_geneventmap() const
  {
    return _geneventmap;
  }

  PHHepMCGenEventMap *get_geneventmap()
  {
    return _geneventmap;
  }

  gsl_rng *get_random_generator()
  {
    return RandomGenerator;
  }

  void set_geneventmap(PHHepMCGenEventMap *geneventmap)
  {
    _geneventmap = geneventmap;
  }

  //! Beam angle in lab polar coordinate.
  //! @param[in] beamA_theta beamA, in pair of Theta-Phi. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamA_phi beamA, in pair of Theta-Phi. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_theta  beamB, in pair of Theta-Phi. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_phi  beamB, in pair of Theta-Phi. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  void set_beam_direction_theta_phi(
      const double beamA_theta,
      const double beamA_phi,
      const double beamB_theta,
      const double beamB_phi)
  {
    m_beam_direction_theta_phi = {{beamA_theta, beamA_phi}, {beamB_theta, beamB_phi}};
  }

  //! Beam angle divergence in accelerator beam coordinate.
  //! @param[in] beamA_divergence_h beamA, horizontal divergence Gaussian Sigma. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamA_divergence_v beamA, vertical divergence Gaussian Sigma. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_divergence_h beamB, horizontal divergence Gaussian Sigma. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_divergence_v beamB, vertical divergence Gaussian Sigma. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  void set_beam_angular_divergence_hv(
      const double beamA_divergence_h,
      const double beamA_divergence_v,
      const double beamB_divergence_h,
      const double beamB_divergence_v)
  {
    m_beam_angular_divergence_hv = {{beamA_divergence_h, beamA_divergence_v}, {beamB_divergence_h, beamB_divergence_v}};
  }

  //! Central beam angle shift as linear function of longitudinal vertex position, d_shift/dz,
  //! which is used to represent leading order effect of crab cavity momentum kick on the beam bunch
  //! @param[in] beamA_h beamA, horizontal angle dh/dz. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamA_v beamA, vertical angle dv/dz. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_h beamB, horizontal angle dh/dz. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_v beamB, vertical angle dv/dz. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  void set_beam_angular_z_coefficient_hv(
      const double beamA_h,
      const double beamA_v,
      const double beamB_h,
      const double beamB_v)
  {
    m_beam_angular_z_coefficient_hv = {{beamA_h, beamA_v}, {beamB_h, beamB_v}};
  }

  //! simulate bunch interaction instead of applying vertex distributions
  void use_beam_bunch_sim(bool b) { m_use_beam_bunch_sim = b; }

  //! Beam bunch geometry as 3D Gauss width
  //! First element is beamA, in vector of Gaussian Sigma H,V,Longitudinal
  //! Second element is beamB, in vector of Gaussian Sigma H,V,Longitudinal
  void set_beam_bunch_width(const std::vector<double> &beamA, const std::vector<double> &beamB);

  void CopySettings(PHHepMCGenHelper &helper);

  //! copy setting to helper_dest
  void CopySettings(PHHepMCGenHelper *helper_dest);

  //! copy setting from helper_src
  void CopyHelperSettings(PHHepMCGenHelper *helper_src);

  void Print(const std::string &what = "ALL") const;

  void PHHepMCGenHelper_Verbosity(int v) { m_verbosity = v; }

  int PHHepMCGenHelper_Verbosity() { return m_verbosity; }

 protected:
  //! Record the translation,boost,rotation for HepMC frame to lab frame according to collision settings
  void HepMC2Lab_boost_rotation_translation(PHHepMCGenEvent *genevent);

  //! move vertex in translation according to vertex settings
  void move_vertex(PHHepMCGenEvent *genevent);

  //! generate vertx with bunch interaction according to
  //! https://github.com/eic/documents/blob/d06b5597a0a89dcad215bab50fe3eefa17a097a5/reports/general/Note-Simulations-BeamEffects.pdf
  //! \return pair of bunch local z position for beam A and beam B
  std::pair<double, double> generate_vertx_with_bunch_interaction(PHHepMCGenEvent *genevent);

 private:
  gsl_rng *RandomGenerator{nullptr};

  double smear(const double position, const double width, VTXFUNC dist) const;

  //! function to convert spherical coordinate to Hep3Vector in x-y-z
  static CLHEP::Hep3Vector pair2Hep3Vector(const std::pair<double, double> &theta_phi);

  VTXFUNC _vertex_func_x{Gaus};
  VTXFUNC _vertex_func_y{Gaus};
  VTXFUNC _vertex_func_z{Gaus};
  VTXFUNC _vertex_func_t{Gaus};

  double _vertex_x{0.};
  double _vertex_y{0.};
  double _vertex_z{0.};
  double _vertex_t{0.};

  double _vertex_width_x{0.};
  double _vertex_width_y{0.};
  double _vertex_width_z{0.};
  double _vertex_width_t{0.};

  //! Beam angle in lab polar coordinate.
  //! First element is beamA, in pair of Theta-Phi. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! Second element is beamB, in pair of Theta-Phi. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  std::pair<std::pair<double, double>, std::pair<double, double>> m_beam_direction_theta_phi = {
      {0, 0},    //+z beam
      {M_PI, 0}  //-z beam
  };

  //! Beam angle divergence in accelerator beam coordinate.
  //! First element is beamA, in pair of Gaussian Sigma_H Sigma_V. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! Second element is beamB, in pair of Gaussian Sigma_H Sigma_V. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  std::pair<std::pair<double, double>, std::pair<double, double>> m_beam_angular_divergence_hv = {
      {0, 0},  //+z beam
      {0, 0}   //-z beam
  };

  //! Central beam angle shift as linear function of longitudinal vertex position, d_shift/dz,
  //! which is used to represent leading order effect of crab cavity momentum kick on the beam bunch
  //! First element is beamA, in pair of dh/dz, dv/dz. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! Second element is beamB, in pair of dh/dz, dv/dz. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  std::pair<std::pair<double, double>, std::pair<double, double>> m_beam_angular_z_coefficient_hv = {
      {0, 0},  //+z beam
      {0, 0}   //-z beam
  };

  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id{0};

  //! whether reuse vertex from another PHHepMCGenEvent
  bool _reuse_vertex{false};

  //! if _reuse_vertex, which embedding_id provide the source vertex. Additional smearing and shift possible with set_vertex_distribution_*()
  int _reuse_vertex_embedding_id{std::numeric_limits<int>::min()};

  //! pointer to the output container
  PHHepMCGenEventMap *_geneventmap{nullptr};

  //! verbosity
  int m_verbosity{0};

  //! simulate bunch interaction instead of applying vertex distributions
  bool m_use_beam_bunch_sim{false};

  //! Beam bunch geometry as 3D Gauss width
  //! First element is beamA, in vector of Gaussian Sigma H,V,Longitudinal
  //! Second element is beamB, in vector of Gaussian Sigma H,V,Longitudinal
  std::pair<std::vector<double>, std::vector<double>> m_beam_bunch_width = {
      {0, 0, 0},  //+z beam
      {0, 0, 0}   //-z beam
  };

  //! use m_beam_bunch_width to calculate horizontal and vertical collision width
  //! \param[in] hv_index 0: horizontal. 1: vertical
  double get_collision_width(unsigned int hv_index);
};

#endif /* PHHEPMC_PHHEPMCGENHELPER_H */
