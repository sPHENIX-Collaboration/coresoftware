#ifndef __PHPYTHIA6_H__
#define __PHPYTHIA6_H__

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenEvent.h>

#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <string>

class PHCompositeNode;
class PHHepMCGenEvent;

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

  void set_config_file( const std::string cfg_file ) { _configFile = cfg_file; }

  void print_config() const;

  void set_node_name(std::string s) {_node_name = s;}

  void save_ascii( std::string fname = "pythia_hepmc.dat" )
  {
    _save_ascii = true;
    _filename_ascii = fname;
  }

private:

  int ReadConfig(const std::string cfg_file = "");
  int CreateNodeTree(PHCompositeNode *topNode);

  int _eventcount;

  // output
  std::string _node_name;

  // Pythia6 configuration file
  std::string _configFile;

  // HepMC
  PHHepMCGenEvent *_phhepmcevt;

  /**
   * Save HepMC event to ASCII file?
   */
  bool _save_ascii;

  /**
   * ASCII file name to save HepMC event to
   */
  std::string _filename_ascii;

  /**
   * definition needed to use pythia wrapper headers from HepMC
   */
  void initPythia();

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif	/* __PHPYTHIA6_H__ */

