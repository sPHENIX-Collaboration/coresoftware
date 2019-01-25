#ifndef PHPYTHIA6_PHPYTHIA6_H
#define PHPYTHIA6_PHPYTHIA6_H

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenHelper.h>

#include <string>
#include <vector>

class PHCompositeNode;
class PHHepMCGenEvent;
class PHPy6GenTrigger; 

namespace HepMC {
  class GenEvent;
};

class PHPythia6: public SubsysReco {

public:

  PHPythia6(const std::string &name = "PHPythia6");
  virtual ~PHPythia6();

  int Init(PHCompositeNode *topNode);

  int process_event(PHCompositeNode *topNode);

  int ResetEvent(PHCompositeNode *topNode);

  int End(PHCompositeNode *topNode);

  void set_config_file( const std::string &cfg_file ) { _configFile = cfg_file; }

  void print_config() const;

  void beam_vertex_parameters(double beamX,
            double beamY,
            double beamZ,
            double beamXsigma,
            double beamYsigma,
            double beamZsigma) {

    set_vertex_distribution_mean(beamX, beamY, beamZ, 0);
    set_vertex_distribution_width(beamXsigma, beamYsigma, beamZsigma, 0);
  }

  void save_ascii( const std::string &fname = "pythia_hepmc.dat" )
  {
    _save_ascii = true;
    _filename_ascii = fname;
  }

  /// set event selection criteria
  void register_trigger(PHPy6GenTrigger *theTrigger);
  void set_trigger_OR() { _triggersOR = true; } // default true
  void set_trigger_AND() { _triggersAND = true; }

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

  int ReadConfig(const std::string &cfg_file = "");
  int CreateNodeTree(PHCompositeNode *topNode);

  /** Certain Pythia switches and parameters only accept integer values
   * This function checks if input values are integers and
   * warns the user if they are not
   */
  void IntegerTest(double number );

  int _eventcount;
  int _geneventcount;

  // Pythia6 configuration file
  std::string _configFile;

  /**
   * Save HepMC event to ASCII file?
   */
  bool _save_ascii;

  /**
   * ASCII file name to save HepMC event to
   */
  std::string _filename_ascii;

  // event selection
  std::vector<PHPy6GenTrigger*> _registeredTriggers;
  bool _triggersOR;
  bool _triggersAND;
 
  /**
   * definition needed to use pythia wrapper headers from HepMC
   */
  void initPythia();
  //! helper for insert HepMC event to DST node and add vertex smearing
  PHHepMCGenHelper hepmc_helper;

};

#endif	/* PHPYTHIA6_PHPYTHIA6_H */

