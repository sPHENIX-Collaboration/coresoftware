#ifndef __PHPYTHIA8_H__
#define __PHPYTHIA8_H__

#include <fun4all/SubsysReco.h>
#include <phhepmc/PHHepMCGenEvent.h>

#ifndef __CINT__
#include <Pythia8/Pythia.h>
#endif

#include <Rtypes.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

class PHCompositeNode;
class PHHepMCGenEvent;
class PHHepMCFilter;
class TTree;
class TFile;
class TRandom;

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
  
  //! constructor
  PHPythia8(const std::string &name = "PHPythia8");
  
  //! destructor
  virtual ~PHPythia8();

  //! configuration file
  void SetConfigFile( const char* cfg_file ) {
    if ( cfg_file ) _configFile = cfg_file;
  }

  //! Read Config File
  /*! if argument is 0 current _configFile is used. It is overwritten otherwise */
  int ReadConfig(const char *cfg_file = 0);

  //! Print Config File
  void PrintConfig() const;

  //! Set Seed of random number generator
  void SetSeed(const int s) { _seed = s; }

  //! Methods Derived from SubsysReco
  int Init(PHCompositeNode *topNode);
  
  //int InitRun(PHCompositeNode *topNode);
  
  //! event method
  int process_event(PHCompositeNode *topNode);
  
  //! event reset
  int ResetEvent(PHCompositeNode *topNode);
  
  //! end of job
  int End(PHCompositeNode *topNode);

  //Triggers
  void registerTrigger(PHPy8GenTrigger *theTrigger);
  void setTriggerOR() { _triggersOR = true; } // default true
  void setTriggerAND() { _triggersAND = true; }

  void ProcessString(std::string s) {_commands.push_back(s);}
  
  void SetNodeName(std::string s) {_node_name = s;}

  void BeamVertexDist(double beamX, double beamY, double beamZ,
		      double beamXsigma, double beamYsigma, double beamZsigma) {
    _useBeamVtx = true;
    _beamX = beamX;
    _beamY = beamY;
    _beamZ = beamZ;    
    _beamXsigma = beamX;
    _beamYsigma = beamY;
    _beamZsigma = beamZ;
  }

private:
  
  int CreateNodeTree(PHCompositeNode *topNode);
  double percentDiff(const double a, const double b){return abs((a-b)/a);}
  
  int _eventcount;

  // output
  std::string _node_name;

  // vertex placement
  TRandom *_rand;
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
  long int _seed;		
  
  // HepMC
  HepMC::Pythia8ToHepMC *_pythiaToHepMC;
  PHHepMCGenEvent *_phhepmcevt;
};

#endif	/* __PHPYTHIA8_H__ */

