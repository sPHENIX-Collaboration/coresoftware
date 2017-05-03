#ifndef __SHEPGEN_H__
#define __SHEPGEN_H__

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenEvent.h>

#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <iostream>
#include <string>

class PHCompositeNode;
class PHHepMCGenEvent;

class HGenManager;
class HParamManager;

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

void set_node_name(std::string s) {_node_name = s;}

private:

  /** Print HEPGen++ logo to screen
   */
  void printlogo();

  /** Create node tree
   */
  int create_node_tree(PHCompositeNode *topNode);

  int _eventcount;

  // output
  std::string _node_name;

  // HEPGen++
  HGenManager* _hgenManager;
  HParamManager* _hgenParManager;

  std::string _datacardFile;

  // HepMC
  PHHepMCGenEvent *_phhepmcevt;

//#ifndef __CINT__
//  gsl_rng *RandomGenerator;
//#endif
};

#endif  /* __SHEPGEN_H__ */

