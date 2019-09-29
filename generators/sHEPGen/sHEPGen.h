#ifndef __SHEPGEN_H__
#define __SHEPGEN_H__

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#if !defined(__CINT__) || defined(__CLING__)
#include <gsl/gsl_rng.h>
#endif

#include <iostream>
#include <string>

class PHCompositeNode;
class PHHepMCGenEvent;

class HGenManager;
class HLorentzVector;

namespace HepMC {
  class GenEvent;
};


class sHEPGen: public SubsysReco {

public:

  sHEPGen(const std::string &name = "sHEPGen");
  virtual ~sHEPGen();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void set_datacard_file( const char* cfg_file ) {
    if ( cfg_file ) _datacardFile = cfg_file;
  }

  void set_momentum_electron( double emom ) {
    _p_electron_lab = emom;
  }

  void set_momentum_hadron( double hmom ) {
    _p_hadron_lab = hmom;
  }

  void beam_vertex_parameters(double beamX,
            double beamY,
            double beamZ,
            double beamXsigma,
            double beamYsigma,
            double beamZsigma) {

    set_vertex_distribution_mean(beamX, beamY, beamZ, 0);
    set_vertex_distribution_width(beamXsigma, beamYsigma, beamZsigma, 0);
  }

  //! toss a new vertex according to a Uniform or Gaus distribution
  void set_vertex_distribution_function(PHHepMCGenHelper::VTXFUNC x, PHHepMCGenHelper::VTXFUNC y, PHHepMCGenHelper::VTXFUNC z, PHHepMCGenHelper::VTXFUNC t)
  {
    hepmc_helper.set_vertex_distribution_function(x, y, z, t);
  }

  //! set the mean value of the vertex distribution, use PHENIX units of cm, ns
  void set_vertex_distribution_mean(const double x, const double y, const double z, const double t)
  {
    hepmc_helper.set_vertex_distribution_mean(x, y, z, t);
  }

  //! set the width of the vertex distribution function about the mean, use PHENIX units of cm, ns
  void set_vertex_distribution_width(const double x, const double y, const double z, const double t)
  {
    hepmc_helper.set_vertex_distribution_width(x, y, z, t);
  }
  //
  //! reuse vertex from another PHHepMCGenEvent with embedding_id = src_embedding_id Additional smearing and shift possible with set_vertex_distribution_*()
  void set_reuse_vertex(int src_embedding_id)
  {
    hepmc_helper.set_reuse_vertex(src_embedding_id);
  }

  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  int get_embedding_id() const { return hepmc_helper.get_embedding_id(); }
  //
  //! embedding ID for the event
  //! positive ID is the embedded event of interest, e.g. jetty event from pythia
  //! negative IDs are backgrounds, .e.g out of time pile up collisions
  //! Usually, ID = 0 means the primary Au+Au collision background
  void set_embedding_id(int id) { hepmc_helper.set_embedding_id(id); }
private:

  /** Print HEPGen++ logo to screen
   */
  void printlogo();

  /** Create node tree
   */
  int create_node_tree(PHCompositeNode *topNode);

  int _eventcount;

  double _p_electron_lab;
  double _p_hadron_lab;

  HLorentzVector *_p4_electron_lab;
  HLorentzVector *_p4_hadron_lab;
  HLorentzVector *_p4_hadron_lab_invert;
  HLorentzVector *_p4_electron_prest;
  HLorentzVector *_p4_hadron_prest;

  // HEPGen++
  HGenManager* _hgenManager;

  std::string _datacardFile;

  //! helper for insert HepMC event to DST node and add vertex smearing
  PHHepMCGenHelper hepmc_helper;
};

#endif  /* __SHEPGEN_H__ */

