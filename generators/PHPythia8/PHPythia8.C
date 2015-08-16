#include "PHPythia8.h"

#include "PHPy8GenTrigger.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>
#include <phool/PHIODataNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHObject.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHNodeReset.h>
#include <phool/PHTimeStamp.h>

#include <TMCParticle.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TRandom.h>

#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC2.h>
#include <HepMC/GenEvent.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdlib.h>
#include <ctime>
#include <sys/time.h>
#include <algorithm>
#include <cctype>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

//static const Float_t CM2MM = 10.; // cm to mm comversion \todo why not use HepMC functions?

PHPythia8::PHPythia8(const std::string &name): 
  SubsysReco(name),
  _eventcount(0),
  _node_name("PHHepMCGenEvent"),
  _rand(NULL),
  _useBeamVtx(false),
  _beamX(0),
  _beamXsigma(0),
  _beamY(0),
  _beamYsigma(0),
  _beamZ(0),
  _beamZsigma(0),
  _registeredTriggers(),
  _triggersOR(true),
  _triggersAND(false),
  _pythia(NULL),
  _configFile("phpythia8.cfg"),
  _commands(),
  _seed(-1),  
  _pythiaToHepMC(NULL),
  _phhepmcevt(NULL) {

  char *charPath = getenv("PYTHIA8");
  if (!charPath) {
    cout << "PHPythia8::Could not find $PYTHIA8 path!" << endl;
    return;
  }
  
  std::string thePath(charPath);
  thePath += "/xmldoc/";
  _pythia = new Pythia8::Pythia(thePath.c_str());

  _pythiaToHepMC = new HepMC::Pythia8ToHepMC();
  _pythiaToHepMC->set_store_proc(true);
  _pythiaToHepMC->set_store_pdf(true);
  _pythiaToHepMC->set_store_xsec(true);  
}

PHPythia8::~PHPythia8() { 
  if (_pythia) delete _pythia;  
  if (_rand) delete _rand;
}

int PHPythia8::Init(PHCompositeNode *topNode) {
  
  if (!_configFile.empty()) read_config();  
  for (unsigned int j = 0; j < _commands.size(); j++) {
    _pythia->readString(_commands[j]);
  }
  
  create_node_tree(topNode);

  // event numbering will start from 1
  _eventcount = 0;

  // PYTHIA8 has very specific requires for its random number range
  // I map the designated unique seed from recoconst into something
  // acceptable for PYTHIA8
  
  if (_seed < 0) {
    _seed = abs(_seed);
  }
  
  while (_seed>900000000) {
    _seed = _seed - 900000000;
  }
  
  if ( (_seed>=0) && (_seed<=900000000) ) {
    _pythia->readString("Random:setSeed = on");
    _pythia->readString(Form("Random:seed = %lu",_seed));
  } else {
    cout << PHWHERE << " ERROR: seed " << _seed << " is not valid" << endl;
    exit(2); 
  }

  _pythia->init();

  print_config();

  if (_useBeamVtx) _rand = new TRandom(_seed);
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::End(PHCompositeNode *topNode) {
  //-* dump out closing info (cross-sections, etc)
  _pythia->stat();
  
  if (verbosity > 1) cout << "PHPythia8::End - I'm here!" << endl;

  //match pythia printout
  cout << " |                                                                "
       << "                                                 | " << endl; 
  cout << "                         PHPythia8::End - " << _eventcount
       << " events passed trigger" << endl;
  cout << "                         Fraction passed: " << _eventcount
       << "/" << _pythia->info.nAccepted()
       << " = " << _eventcount/float(_pythia->info.nAccepted()) << endl;
  cout << " *-------  End PYTHIA Trigger Statistics  ------------------------"
       << "-------------------------------------------------* " << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________________________________________
int PHPythia8::read_config(const char *cfg_file) {

  if ( cfg_file ) _configFile = cfg_file;
  cout << "PHPythia8::read_config - Reading " << _configFile << endl;
  
  ifstream infile( _configFile.c_str() ); 
  if (infile.fail ()) {
    cout << "PHPythia8::read_config - Failed to open file " << _configFile << endl;    
    exit(2);
  }

  _pythia->readFile(_configFile.c_str());

  return Fun4AllReturnCodes::EVENT_OK;
}

//-* print pythia config info
void PHPythia8::print_config() const {
  _pythia->info.list();
  cout << "Using seed " << _seed << endl;
}

int PHPythia8::process_event(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "PHPythia8::process_event - event: " << _eventcount << endl;
  
  bool passedGen = false;
  bool passedTrigger = false;
  std::vector<bool> theTriggerResults;
  int genCounter = 0;

  while (!passedTrigger) {
    ++genCounter;

    // generate another pythia event
    while (!passedGen) {
      passedGen = _pythia->next();
    }

    // test trigger logic
    
    bool andScoreKeeper = true;
    if (verbosity > 2) {
      cout << "PHPythia8::process_event - triggersize: " << _registeredTriggers.size() << endl;
    }

    for (unsigned int tr = 0; tr < _registeredTriggers.size(); tr++) { 
      bool trigResult = _registeredTriggers[tr]->Apply(_pythia);

      if (verbosity > 2) {
	cout << "PHPythia8::process_event trigger: "
	     << _registeredTriggers[tr]->GetName() << "  " << trigResult << endl;
      }

      if (_triggersOR && trigResult) {
	passedTrigger = true;
	break;
      } else if (_triggersAND) {
	andScoreKeeper &= trigResult;
      }
      
      if (verbosity > 2 && !passedTrigger) {
	cout << "PHPythia8::process_event - failed trigger: "
	     << _registeredTriggers[tr]->GetName() <<  endl;
      }
    }

    if ((andScoreKeeper && _triggersAND) || (_registeredTriggers.size() == 0)) {
      passedTrigger = true;
      genCounter = 0;
    }

    passedGen = false;
  }

  // fill HepMC object with event & pass to 
  
  HepMC::GenEvent *genevent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
  _pythiaToHepMC->fill_next_event(*_pythia, genevent, _eventcount);

  // pass HepMC to PHNode
  
  bool success = _phhepmcevt->addEvent(genevent);
  if (!success) {
    cout << "PHPythia8::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // shift node if needed  
  if (_useBeamVtx) {
    double mvVtxX = _rand->Gaus(_beamX,_beamXsigma);
    double mvVtxY = _rand->Gaus(_beamY,_beamYsigma);
    double mvVtxZ = _rand->Gaus(_beamZ,_beamZsigma);
    _phhepmcevt->moveVertex(mvVtxX,mvVtxY,mvVtxZ,0.0);
  }

  // print outs
  
  if (verbosity > 2) cout << "PHPythia8::process_event - FINISHED WHOLE EVENT" << endl;
  if (_eventcount < 2 && verbosity > 1) _pythia->event.list();
  if (_eventcount >= 2 && verbosity > 5) _pythia->event.list();

  ++_eventcount;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::create_node_tree(PHCompositeNode *topNode) {

  PHCompositeNode *dstNode;
  PHNodeIterator iter(topNode);

  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing doing nothing" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _phhepmcevt = new PHHepMCGenEvent();
  PHObjectNode_t *newNode = new PHObjectNode_t(_phhepmcevt,_node_name.c_str(),"PHObject");
  dstNode->addNode(newNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::ResetEvent(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHPythia8::register_trigger(PHPy8GenTrigger *theTrigger) {
  if(verbosity > 1) cout << "PHPythia8::registerTrigger - trigger " << theTrigger->GetName() << " registered" << endl;
  _registeredTriggers.push_back(theTrigger);
}
