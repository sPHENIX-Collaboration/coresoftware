// $Id: $

/*!
 * \file PHHepMCGenHelper.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_PHHEPMC_PHHEPMCGENHELPER_H_
#define SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_PHHEPMC_PHHEPMCGENHELPER_H_



class PHCompositeNode;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;

namespace HepMC
{
class GenEvent;
};

#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

/*!
 * \brief PHHepMCGenHelper
 */
class PHHepMCGenHelper
{
 public:
  PHHepMCGenHelper();
  virtual ~PHHepMCGenHelper();

  //! supported function distributions
  enum VTXFUNC
  {
    Uniform,
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
  //! init interface nodes
  int create_node_tree(PHCompositeNode *topNode);

  //! send HepMC::GenEvent to DST tree. This function takes ownership of evt
  PHHepMCGenEvent *insert_event(HepMC::GenEvent *evt);

  //! move vertex according to vertex settings
  void move_vertex(PHHepMCGenEvent *genevent);

  const PHHepMCGenEventMap *get_geneventmap() const
  {
    return _geneventmap;
  }

  PHHepMCGenEventMap *get_geneventmap()
  {
    return _geneventmap;
  }

  void set_geneventmap(PHHepMCGenEventMap *geneventmap)
  {
    _geneventmap = geneventmap;
  }

 protected:
#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif

  double smear(const double position, const double width, VTXFUNC dist) const;

 private:
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

  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int _embedding_id;

  PHHepMCGenEventMap *_geneventmap;
};

#endif /* SIMULATION_CORESOFTWARE_SIMULATION_G4SIMULATION_PHHEPMC_PHHEPMCGENHELPER_H_ */
