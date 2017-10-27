#include "PHPythia8.h"

#include "PHPy8GenTrigger.h"

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHIODataNode.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>


#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC2.h>
#include <HepMC/GenEvent.h>

#include <gsl/gsl_randist.h>

#include <TString.h> // needed for Form()

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

//static const Float_t CM2MM = 10.; // cm to mm comversion \todo why not use HepMC functions?

PHPythia8::PHPythia8(const std::string &name): 
  SubsysReco(name),
  _eventcount(0),
  _registeredTriggers(),
  _triggersOR(true),
  _triggersAND(false),
  _pythia(NULL),
  _configFile("phpythia8.cfg"),
  _commands(),
  _pythiaToHepMC(NULL) {

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

  hepmc_helper.set_embedding_id(1); // default embedding ID to 1
}

PHPythia8::~PHPythia8() {
  delete _pythia;  
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

  unsigned int seed = PHRandomSeed();

  if (seed > 900000000) {
    seed = seed % 900000000;
  }

  if ( (seed>0) && (seed<=900000000) ) {
    _pythia->readString("Random:setSeed = on");
    _pythia->readString(Form("Random:seed = %u",seed));
  } else {
    cout << PHWHERE << " ERROR: seed " << seed << " is not valid" << endl;
    exit(1); 
  }

  _pythia->init();

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
  //_pythia->info.list();
}

int PHPythia8::process_event(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "PHPythia8::process_event - event: " << _eventcount << endl;
  
  bool passedGen = false;
  bool passedTrigger = false;
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

  
  /* pass HepMC to PHNode*/
  PHHepMCGenEvent * success = hepmc_helper . insert_event(genevent);
  if (!success) {
    cout << "PHPythia8::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // print outs
  
  if (verbosity > 2) cout << "PHPythia8::process_event - FINISHED WHOLE EVENT" << endl;
  if (_eventcount < 2 && verbosity > 1) _pythia->event.list();
  if (_eventcount >= 2 && verbosity > 5) _pythia->event.list();

  ++_eventcount;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::create_node_tree(PHCompositeNode *topNode) {

  hepmc_helper.create_node_tree(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::ResetEvent(PHCompositeNode *topNode) {
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHPythia8::register_trigger(PHPy8GenTrigger *theTrigger) {
  if(verbosity > 1) cout << "PHPythia8::registerTrigger - trigger " << theTrigger->GetName() << " registered" << endl;
  _registeredTriggers.push_back(theTrigger);
}
