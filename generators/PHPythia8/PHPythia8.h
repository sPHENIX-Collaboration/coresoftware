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
  void SetSeed(const int s) { fSeed = s; }

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

  void GaussianVertexDist(const double mean = 0.0, const double sigma = 15.0) {
    _useGaussianVtx = true;
    _gaussMean = mean;
    _gaussSigma = sigma;
  }

  void BeamDiamondVtx(int Run = 12) {
    _useBeamVtx = true;
    if (Run == 12) {
      _beamX = 0.1247;
      _beamY = -0.1601;
      _beamZ = 2.478;

      _beamXsigma = .021;
      _beamYsigma = .018;
      _beamZsigma = 15.75;
    } else {
      std::cout << "PHPythia8::BeamDiamondVtx - Unknown run " << Run << ". Using Run 12..." << std::endl;
      BeamDiamondVtx(12);
    }
  }

  void CorrelateQ(std::string s) {
    std::cout << std::endl;
    std::cout << "PHPythia8::CorrelateQ - Will correlate Q with event on node " << s << "!!!" << std::endl;
    std::cout << std::endl;
    _qNodeName = s; 
    _correlateQ = true;
  }

  double percentDiff(const double a, const double b){return abs((a-b)/a);}

private:
  
  //! node tree
  int CreateNodeTree(PHCompositeNode *topNode);
  
  //! event
  int eventcount;
  
  //! configuration file. Default is "phpythia.cfg
  std::string _configFile;
  
  //! pythia interface
  #ifndef __CINT__
  Pythia8::Pythia *pythia;
  #endif
  
  //! seed to random number generator
  long int fSeed;		
  
  //HepMC
  HepMC::GenEvent *hepmcevt;
  HepMC::Pythia8ToHepMC *pythiaToHepMC;
  bool _isHepMC;
  bool _isBoth;
  bool _isPHPythia;
  PHHepMCGenEvent *phhepmcevt;

  std::vector<PHPy8GenTrigger*> _registeredTriggers;
  std::vector<std::string> _commands;

  bool _triggersOR;
  bool _triggersAND;
  
  std::string _node_name;
  
  bool _useGaussianVtx;
  double _gaussMean, _gaussSigma;
  TRandom *rand;

  std::string _qNodeName;
  bool _correlateQ;

  bool _useBeamVtx;
  double _beamX, _beamXsigma;
  double _beamY, _beamYsigma;
  double _beamZ, _beamZsigma;

};

#endif	/* __PHPYTHIA8_H__ */

