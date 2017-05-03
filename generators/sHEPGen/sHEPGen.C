#include "sHEPGen.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHIODataNode.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>

/* HEPGen includes */
#include <hgenmanager.h>
//#include "hvector.h"
//#include "config.h"
//#include "hlorentzvector.h"
//#include "hcardparser.h"
//#include "hpionicdata.h"

/* HepMC includes */
#include <HepMC/GenEvent.h>

/* ROOT includes */
#include <TString.h> // needed for Form()

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;


sHEPGen::sHEPGen(const std::string &name):
  SubsysReco(name),
  _eventcount(0),
  _p_electron_lab(20),
  _p_hadron_lab(250),
  _detailed_debug(false),
  _node_name("PHHepMCGenEvent"),
  _hgenManager(NULL),
  _datacardFile("hepgen_dvcs.data"),
  _phhepmcevt(NULL)
{

}

sHEPGen::~sHEPGen() {

}

int sHEPGen::Init(PHCompositeNode *topNode) {

  printlogo();

  /* get instance of HepGenManager */
  _hgenManager = HGenManager::getInstance();
  _hgenManager->loadSettings( _datacardFile );
  _hgenManager->setupGenerator();

  /* set beam parameters */
  double beamE = 165.47143;
  cout << "----ELEPT: " << beamE << endl;
  _hgenManager->getParamManager()->getStruct()->ELEPT = beamE;
  _hgenManager->getParamManager()->getStruct()->PARL.at(2) = beamE;

  /* set random seed */
  unsigned int seed = PHRandomSeed();
  _hgenManager->setSeed ( seed );

  /* enable detailed event record printput for debugging */
  if ( _detailed_debug )
    _hgenManager->enableDebug();

  create_node_tree(topNode);

  _eventcount = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int sHEPGen::End(PHCompositeNode *topNode) {

  cout << "Reached the sHEPGen::End()" << endl;

  //-* dump out closing info (cross-sections, etc)

  if (verbosity > 1) cout << "sHEPGen::End - I'm here!" << endl;



  return Fun4AllReturnCodes::EVENT_OK;
}

//__________________________________________________________
//int sHEPGen::read_config(const char *cfg_file) {
//
//  if ( cfg_file ) _configFile = cfg_file;
//  cout << "sHEPGen::read_config - Reading " << _configFile << endl;
//
//  ifstream infile( _configFile.c_str() );
//  if (infile.fail ()) {
//    cout << "sHEPGen::read_config - Failed to open file " << _configFile << endl;
//    exit(2);
//  }
//
//  _pythia->readFile(_configFile.c_str());
//
//  return Fun4AllReturnCodes::EVENT_OK;
//}
//
////-* print pythia config info
//void sHEPGen::print_config() const {
//  //_pythia->info.list();
//}

int sHEPGen::process_event(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "sHEPGen::process_event - event: " << _eventcount << endl;

  _hgenManager->oneShot();
//  _hgenManager->getEvent()->printDebug();


  // fill HepMC object with event & pass to
  //  HepMC::GenEvent *genevent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);

  // pass HepMC to PHNode
  //  bool success = _phhepmcevt->addEvent(genevent);
  //  if (!success) {
  //    cout << "sHEPGen::process_event - Failed to add event to HepMC record!" << endl;
  //    return Fun4AllReturnCodes::ABORTRUN;
  //  }

  // print outs
  if (verbosity > 2) cout << "sHEPGen::process_event - FINISHED WHOLE EVENT" << endl;

  ++_eventcount;

  return Fun4AllReturnCodes::EVENT_OK;
}


int sHEPGen::create_node_tree(PHCompositeNode *topNode) {

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


void sHEPGen::printlogo() {

  cout << endl <<  endl;
  for ( int i = 0; i < hepconst::logolength; i++ )
    cout << hepconst::logo[i] << endl;
  cout << endl <<  endl;

}
