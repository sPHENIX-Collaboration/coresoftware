#include "sHEPGen.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHIODataNode.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>

#include <gsl/gsl_randist.h>

/* HEPGen includes */
#include <hgenmanager.h>
#include <hparammanager.h>
//#include "hvector.h"
//#include "config.h"
//#include "hlorentzvector.h"
//#include "hcardparser.h"
//#include "hpionicdata.h"
//#include "hbooker.h"
//#include "hbookbackendASCII.h"
//#ifdef USE_ROOT
//#include "hbookbackendROOT.h"
//#endif

/* HepMC includes */
#include <HepMC/GenEvent.h>

/* ROOT includes */
#include <TString.h> // needed for Form()

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;


sHEPGen::sHEPGen(const std::string &name):
  SubsysReco(name),
  _eventcount(0),
  _node_name("PHHepMCGenEvent"),
  _hgenManager(NULL),
  _hgenParManager(NULL),
  _datacardFile("hepgen_dvcs.data"),
  _phhepmcevt(NULL) {

  //  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
}

sHEPGen::~sHEPGen() {
  //  gsl_rng_free (RandomGenerator);
  delete _hgenParManager;
}

int sHEPGen::Init(PHCompositeNode *topNode) {

  printlogo();

  cout << "Create manager:" << endl;
  _hgenManager = HGenManager::getInstance();
  cout << _hgenManager << endl;

  cout << "Create parameter manager:" << endl;
  _hgenParManager = new HParamManager( _datacardFile );
  cout << _hgenParManager << endl;

  create_node_tree(topNode);

  _eventcount = 0;

  /* set random seed */
  //  unsigned int seed = PHRandomSeed();
  //  _hgenManager->setSeed ( seed );

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

  double beamE = 165.47143;
  cout << "----ELEPT: " << beamE << endl;
//  ParamManager.getStruct()->ELEPT = beamE;
//  ParamManager.getStruct()->PARL.at(2) = beamE;
//
//  tmp.oneShot();
//
//  tmp.getEvent()->printDebug();


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
