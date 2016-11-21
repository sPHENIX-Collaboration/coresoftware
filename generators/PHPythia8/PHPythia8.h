#ifndef __PHPYTHIA8_H__
#define __PHPYTHIA8_H__

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenEvent.h>

#ifndef __CINT__
#include <Pythia8/Pythia.h>
#endif

#ifndef __CINT__
#include <gsl/gsl_rng.h>
#endif

#include <iostream>
#include <string>

class PHCompositeNode;
class PHHepMCGenEvent;
class PHHepMCFilter;

class PHPy8GenTrigger;

namespace HepMC {
  class GenEvent;
  class Pythia8ToHepMC;
};

namespace Pythia8 {
  class Pythia;
};

class PHPythia8: public SubsysReco {
  
public:
  
  PHPythia8(const std::string &name = "PHPythia8");
  virtual ~PHPythia8();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode); 
  int ResetEvent(PHCompositeNode *topNode); 
  int End(PHCompositeNode *topNode);
  
  void set_config_file( const char* cfg_file ) {
    if ( cfg_file ) _configFile = cfg_file;
  }

  void print_config() const;

  /// set event selection criteria
  void register_trigger(PHPy8GenTrigger *theTrigger);
  void set_trigger_OR() { _triggersOR = true; } // default true
  void set_trigger_AND() { _triggersAND = true; }

  /// pass commands directly to PYTHIA8
  void process_string(std::string s) {_commands.push_back(s);}
  
  void set_node_name(std::string s) {_node_name = s;}

  void beam_vertex_parameters(double beamX,
			      double beamY,
			      double beamZ,
			      double beamXsigma,
			      double beamYsigma,
			      double beamZsigma) {
    _useBeamVtx = true;
    _beamX = beamX;
    _beamY = beamY;
    _beamZ = beamZ;    
    _beamXsigma = beamXsigma;
    _beamYsigma = beamYsigma;
    _beamZsigma = beamZsigma;
  }

private:

  int read_config(const char *cfg_file = 0);
  int create_node_tree(PHCompositeNode *topNode);
  double percent_diff(const double a, const double b){return abs((a-b)/a);}
  
  int _eventcount;

  // output
  std::string _node_name;

  // vertex placement
  bool _useBeamVtx;
  double _beamX, _beamXsigma;
  double _beamY, _beamYsigma;
  double _beamZ, _beamZsigma;

  // event selection
  std::vector<PHPy8GenTrigger*> _registeredTriggers;
  bool _triggersOR;
  bool _triggersAND;
  
  // PYTHIA  
  #ifndef __CINT__
  Pythia8::Pythia *_pythia;
  #endif

  std::string _configFile;
  std::vector<std::string> _commands;
  
  // HepMC
  HepMC::Pythia8ToHepMC *_pythiaToHepMC;
  PHHepMCGenEvent *_phhepmcevt;

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif	/* __PHPYTHIA8_H__ */

