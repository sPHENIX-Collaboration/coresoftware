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

#include <gsl/gsl_rng.h>

#include <cmath>
#include <string>
#include <utility>

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
  const PHHepMCGenEvent * get_PHHepMCGenEvent_template() const;

  //! send HepMC::GenEvent to DST tree. This function takes ownership of evt
  PHHepMCGenEvent *insert_event(HepMC::GenEvent *evt);

  //! Record the translation,boost,rotation for HepMC frame to lab frame according to collision settings
  void HepMC2Lab_boost_rotation_translation(PHHepMCGenEvent *genevent);

  //! move vertex in translation according to vertex settings
  void move_vertex(PHHepMCGenEvent *genevent);

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
  //! @param[in] beamA_divergence_x beamA, in pair of Gaussian Sigma_X Sigma_Y. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamA_divergence_y beamA, in pair of Gaussian Sigma_X Sigma_Y. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_divergence_x beamB, in pair of Gaussian Sigma_X Sigma_Y. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  //! @param[in] beamB_divergence_y beamB, in pair of Gaussian Sigma_X Sigma_Y. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  void set_beam_angular_divergence_xy(
      const double beamA_divergence_x,
      const double beamA_divergence_y,
      const double beamB_divergence_x,
      const double beamB_divergence_y)
  {
    m_beam_angular_divergence_xy = {{beamA_divergence_x, beamA_divergence_y}, {beamB_divergence_x, beamB_divergence_y}};
  }

  void CopySettings(PHHepMCGenHelper &helper);

  //! copy setting to helper_dest
  void CopySettings(PHHepMCGenHelper * helper_dest) ;

  //! copy setting from helper_src
  void CopyHelperSettings(PHHepMCGenHelper * helper_src) ;

  void Print(const std::string &what = "ALL") const;

 private:
  gsl_rng *RandomGenerator;

  double smear(const double position, const double width, VTXFUNC dist) const;

  VTXFUNC _vertex_func_x;
  VTXFUNC _vertex_func_y;
  VTXFUNC _vertex_func_z;
  VTXFUNC _vertex_func_t;

  double _vertex_x;
  double _vertex_y;
  double _vertex_z;
  double _vertex_t;

  double _vertex_width_x;
  double _vertex_width_y;
  double _vertex_width_z;
  double _vertex_width_t;

  //! Beam angle in lab polar coordinate.
  //! First element is beamA, in pair of Theta-Phi. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! Second element is beamB, in pair of Theta-Phi. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  std::pair<std::pair<double, double>, std::pair<double, double>> m_beam_direction_theta_phi = {
      {0, 0},    //+z beam
      {M_PI, 0}  //-z beam
  };

  //! Beam angle divergence in accelerator beam coordinate.
  //! First element is beamA, in pair of Gaussian Sigma_X Sigma_Y. BeamA is aimed to +z direction in the HepMC event generator's coordinate
  //! Second element is beamB, in pair of Gaussian Sigma_X Sigma_Y. BeamA is aimed to -z direction in the HepMC event generator's coordinate
  std::pair<std::pair<double, double>, std::pair<double, double>> m_beam_angular_divergence_xy = {
      {0, 0},  //+z beam
      {0, 0}   //-z beam
  };

  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id;

  //! whether reuse vertex from another PHHepMCGenEvent
  bool _reuse_vertex;

  //! if _reuse_vertex, which embedding_id provide the source vertex. Additional smearing and shift possible with set_vertex_distribution_*()
  int _reuse_vertex_embedding_id;

  //! pointer to the output container
  PHHepMCGenEventMap *_geneventmap;
};

#endif /* PHHEPMC_PHHEPMCGENHELPER_H */
