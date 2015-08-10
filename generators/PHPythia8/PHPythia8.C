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

static const Float_t CM2MM = 10.; // cm to mm comversion \todo why not use HepMC functions?

PHPythia8::PHPythia8(const std::string &name): 
  SubsysReco(name),
  eventcount( 0 ),
  _configFile( "" ),
  fSeed(-1),
  _node_name("PHHepMCGenEvent"),
  _useGaussianVtx(false),
  _gaussMean(0.0),
  _gaussSigma(15.0),
  rand(NULL),
  _qNodeName(""),
  _correlateQ(false),
  _useBeamVtx(false),
  _beamX(0),
  _beamXsigma(0),
  _beamY(0),
  _beamYsigma(0),
  _beamZ(0),
  _beamZsigma(0) {
  
  std::string thePath = getenv("PYTHIA8");
  if (thePath==NULL) cout << "Could not find $PYTHIA8 path!" << endl;
  thePath += "/xmldoc/";
  pythia = new Pythia8::Pythia(thePath.c_str());

  pythiaToHepMC = new HepMC::Pythia8ToHepMC();
  pythiaToHepMC->set_store_proc(true);
  pythiaToHepMC->set_store_pdf(true);
  pythiaToHepMC->set_store_xsec(true);
  
  hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::CM);

  _triggersOR = true;
  _triggersAND = false;
}

PHPythia8::~PHPythia8() { 
  if (pythia!=0) delete pythia;  
  if (_useGaussianVtx) delete rand;
}

int PHPythia8::Init(PHCompositeNode *topNode) {
  
  if (!_configFile.empty()) ReadConfig();  
  for (unsigned int j = 0; j < _commands.size(); j++) pythia->readString(_commands[j]);

  CreateNodeTree(topNode);

  // event numbering will start from 1
  eventcount = 0;

  /// \todo replace with RAND seed from recoconst
  if ( fSeed < 0 ) {
    // first try getting seed from /dev/random
    ifstream devrandom;
    devrandom.open("/dev/random",ios::binary);
    devrandom.read((char*)&fSeed,sizeof(fSeed));
    devrandom.close();
    
    if ( fSeed != -1 ) {
      cout << PHWHERE << " Got seed from /dev/random" << endl;
      fSeed = abs(fSeed)%900000000;
    } else {
      // /dev/random failed, get the random seed from the time of day, to the microsecond
      //fSeed = (Int_t)(time(NULL)/3);
      cout << PHWHERE << " Getting seed from gettimeofday()" << endl;
      timeval xtime;
      int status = gettimeofday(&xtime,NULL);
      if ( status==0 ) {
        fSeed = ((xtime.tv_sec << 12) + (xtime.tv_usec&0xfff))%900000000;
      } else {
        cout << PHWHERE << " something wrong with gettimeofday()" << endl;
      }
    }
  }
  
  if ( (fSeed>=0) && (fSeed<=900000000) ) {
    pythia->readString("Random:setSeed = on");
    pythia->readString(Form("Random:seed = %lu",fSeed));
  } else {
    cout << PHWHERE << " ERROR: seed " << fSeed << " is not valid" << endl;
    exit(2); 
  }

  pythia->init();

  PrintConfig();

  if (_useGaussianVtx || _useBeamVtx) rand = new TRandom(fSeed);
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::End(PHCompositeNode *topNode) {
  //-* dump out closing info (cross-sections, etc)
  pythia->stat();
  
  if (verbosity > 1) cout << "PHPythia8::End - I'm here!" << endl;

  cout << " |                                                                                                                 | " << endl; //match pythia printout
  cout << "                         PHPythia8::End - " << eventcount << " events passed trigger" << endl;
  cout << "                         Fraction passed: " << eventcount << "/" << pythia->info.nAccepted() <<" = " << eventcount/float(pythia->info.nAccepted()) << endl;
  cout << " *-------  End PYTHIA Trigger Statistics  -------------------------------------------------------------------------* " << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________________________________________
int PHPythia8::ReadConfig(const char *cfg_file) {

  if ( cfg_file ) _configFile = cfg_file;
  cout << "PHPythia8::ReadConfig - Reading " << _configFile << endl;
  
  ifstream infile( _configFile.c_str() ); 
  if (infile.fail ()) {
    cout << "PHPythia8::ReadConfig - Failed to open file " << _configFile << endl;
    exit(2);
  }

  pythia->readFile(_configFile.c_str());

  return Fun4AllReturnCodes::EVENT_OK;
}

//-* print pythia config info
void PHPythia8::PrintConfig() const {
  //pythia->init.showProcesses();
  pythia->info.list();
  cout << "Using seed " << fSeed << endl;
}

int PHPythia8::process_event(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "PHPythia8::process_event - event: " << eventcount << endl;

  double _lastEventQ = 0;
  double _lastEventVtxX = 0, _lastEventVtxY = 0, _lastEventVtxZ = 0;
  if (_correlateQ) {
    PHHepMCGenEvent* QhepmcEvent = findNode::getClass<PHHepMCGenEvent>(topNode,_qNodeName.c_str());
    if (!QhepmcEvent) {
      cout << "PHHepMCFilter::process_event - unable to get PHHepMCGenEvent named "
	   << _node_name << ", is Node missing?" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    HepMC::GenEvent *Qevent = QhepmcEvent->getEvent();
    _lastEventQ = Qevent->pdf_info()->scalePDF();
    HepMC::GenEvent::vertex_const_iterator pv = Qevent->vertices_begin();
    _lastEventVtxX =  (*pv)->position().x();
    _lastEventVtxY =  (*pv)->position().y();
    _lastEventVtxZ =  (*pv)->position().z();
  }
  
  bool passedGen = false;
  bool passedTrigger = false;
  std::vector<bool> theTriggerResults;
  int genCounter = 0;
  while (!passedTrigger) {
    genCounter++;
    while (!passedGen) {
      passedGen = pythia->next();
    }
    bool andScoreKeeper = true;
    if (verbosity > 2) cout << "PHPythia8::process_event - triggersize: " << _registeredTriggers.size() << endl;
    for (unsigned int tr = 0; tr < _registeredTriggers.size(); tr++)
	{ 
	  bool trigResult = _registeredTriggers[tr]->Apply(pythia);
	  if(verbosity > 2) cout << "PHPythia8::process_event trigger: " << _registeredTriggers[tr]->GetName() << "  " << trigResult << endl;
	  if(_triggersOR && trigResult)
	    {
	      passedTrigger = true;
	      break;
	    }
	  else if(_triggersAND)
	    {
	      andScoreKeeper &= trigResult;
	    }
	  if(verbosity > 2 && !passedTrigger) cout << "PHPythia8::process_event - failed trigger: " << _registeredTriggers[tr]->GetName() <<  endl;
	}

      if(_correlateQ)
	{
	  double pyQ = pythia->info.QRen();
	  double tolerance = 0.30;
	  if(_lastEventQ < 2)
	    {
	      if(genCounter > 100)   tolerance = 0.45;
	      if(genCounter > 1000)  tolerance = 0.55;
	      if(genCounter > 10000) tolerance = 0.75;
	      if(genCounter > 100000) tolerance = 0.95;
	    }
	  else{
	    if(genCounter > 100)   tolerance = 0.4;
	    if(genCounter > 1000)  tolerance = 0.5;
	    if(genCounter > 10000) tolerance = 0.6;
	    if(genCounter > 100000) tolerance = 0.85;
	  }

	  //std::cout << "Q's: " << _lastEventQ << "   " << pyQ << "    " << percentDiff(_lastEventQ,pyQ) << "  " << tolerance << "  " << genCounter << std::endl;
	  if( percentDiff(_lastEventQ,pyQ) < tolerance && pyQ > 0 )
	    {
	      passedTrigger = true;
	      //std::cout << "Found Q: " << _lastEventQ << "   " << pyQ << "    " << percentDiff(_lastEventQ,pyQ) << "  " << tolerance << "  " << genCounter << std::endl; 
	      genCounter = 0;
	    }
	}  
      if((andScoreKeeper && _triggersAND) || (_registeredTriggers.size() == 0 && !_correlateQ))
	{
	  passedTrigger = true;
	  genCounter = 0;
	}

      passedGen = false;
    }

  if (hepmcevt) delete hepmcevt;
  hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::CM);
  pythiaToHepMC->fill_next_event(*pythia, hepmcevt, eventcount);
  if (!phhepmcevt->addEvent(hepmcevt)) cout << "PHPythia8::process_event - Failed to add event to HepMC record!" << endl;
  if (_useGaussianVtx && !_correlateQ) {
    double mvVtxZ = rand->Gaus(_gaussMean,_gaussSigma);
    phhepmcevt->moveVertex(0.0,0.0,mvVtxZ,0.0);
  }
  if (_useBeamVtx && !_correlateQ) {
    double mvVtxX = rand->Gaus(_beamX,_beamXsigma);
    double mvVtxY = rand->Gaus(_beamY,_beamYsigma);
    double mvVtxZ = 100;
    while (abs(mvVtxZ) > 20) {
      mvVtxZ = rand->Gaus(_beamZ,_beamZsigma);
    }
    phhepmcevt->moveVertex(mvVtxX,mvVtxY,mvVtxZ,0.0);
  }
  if (_correlateQ) {
    phhepmcevt->moveVertex(_lastEventVtxX,_lastEventVtxY,_lastEventVtxZ,0.0);
  }
  
  if (verbosity > 2) cout << "PHPythia8::process_event - FINISHED WHOLE EVENT" << endl;
  if (eventcount < 2 && verbosity > 1) pythia->event.list();   // list full pythia generated event
  if (eventcount >= 2 && verbosity > 5) pythia->event.list();

  ++eventcount;
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHPythia8::CreateNodeTree(PHCompositeNode *topNode) {

  PHCompositeNode *dstNode;
  PHNodeIterator iter(topNode);

  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing doing nothing" << endl;
    return -1;
  }

  phhepmcevt = new PHHepMCGenEvent();
  PHObjectNode_t *newNode = new PHObjectNode_t(phhepmcevt,_node_name.c_str(),"PHObject");
  dstNode->addNode(newNode);

  return 0;
}

int PHPythia8::ResetEvent(PHCompositeNode *topNode) {
  return 0;
}

void PHPythia8::registerTrigger(PHPy8GenTrigger *theTrigger) {
  if(verbosity > 1) cout << "PHPythia8::registerTrigger - trigger " << theTrigger->GetName() << " registered" << endl;
  _registeredTriggers.push_back(theTrigger);
}
