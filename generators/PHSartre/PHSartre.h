#ifndef __PHSARTRE_H__
#define __PHSARTRE_H__

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

class Sartre; 
class Event; 
class EventGeneratorSettings; 
class PHSartreGenTrigger;
class TGenPhaseSpace; 

namespace HepMC {
  class GenEvent;
};

class PHSartre: public SubsysReco {
  
public:
  
  PHSartre(const std::string &name = "PHSartre");
  virtual ~PHSartre();

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode); 
  int ResetEvent(PHCompositeNode *topNode); 
  int End(PHCompositeNode *topNode);
  
  void set_config_file( const char* cfg_file ) {
    if ( cfg_file ) _configFile = cfg_file;
  }

  void print_config() const;

  /// set event selection criteria
  void register_trigger(PHSartreGenTrigger *theTrigger);
  void set_trigger_OR() { _triggersOR = true; _triggersAND = false; } // default true
  void set_trigger_AND() { _triggersAND = true; _triggersOR = false; }

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

  int create_node_tree(PHCompositeNode *topNode);
  double percent_diff(const double a, const double b){return abs((a-b)/a);}
  void randomlyReverseBeams(Event* myEvent);
  void ReverseBeams(Event* myEvent);
  
  int _eventcount;
  int _gencount; 

  // output
  std::string _node_name;

  // vertex placement
  bool _useBeamVtx;
  double _beamX, _beamXsigma;
  double _beamY, _beamYsigma;
  double _beamZ, _beamZsigma;

  // event selection
  std::vector<PHSartreGenTrigger*> _registeredTriggers;
  bool _triggersOR;
  bool _triggersAND;
  
  std::string _configFile;
  std::vector<std::string> _commands;
  
  // HepMC
  PHHepMCGenEvent *_phhepmcevt;

  // Sartre
  Sartre *_sartre;
  EventGeneratorSettings* settings;   
  TGenPhaseSpace *decay; 
  int daughterID;
  double daughterMasses[2];
  bool doPerformDecay;
  

#ifndef __CINT__
  gsl_rng *RandomGenerator;
#endif
};

#endif	/* __PHSARTRE_H__ */

