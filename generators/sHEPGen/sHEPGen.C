#include "sHEPGen.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHIODataNode.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>

/* HEPGen includes */
#include <hgenmanager.h>

/* HepMC includes */
#include <HepMC/GenEvent.h>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;


sHEPGen::sHEPGen(const std::string &name):
  SubsysReco(name),
  _eventcount(0),
  _p_electron_lab(-20),
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

  /* electron and proton mass */
  double mass_e = 5.109989e-4;
  double mass_p = 9.382720e-1;

  /* 4-Vectors of colliding electron and proton in laboratory frame */
  HLorentzVector *p4_electron_lab = new HLorentzVector( 0.,
                                                        0.,
                                                        _p_electron_lab,
                                                        sqrt(_p_electron_lab*_p_electron_lab+mass_e*mass_e) );

  HLorentzVector *p4_hadron_lab = new HLorentzVector( 0.,
                                                      0.,
                                                      _p_hadron_lab,
                                                      sqrt(_p_hadron_lab*_p_hadron_lab+mass_p*mass_p) );

  /* The current version of HEPGen supports only a "Fixed Target" mode, i.e.
     the target (proton) is assumed at rest. Therefore, need to boost collision
     from the laboratory frame to the proton-at-rest frame. */
  HLorentzVector *p4_hadron_lab_invert = new HLorentzVector( 0.,
                                                             0.,
                                                             -1 * _p_hadron_lab,
                                                             sqrt(_p_hadron_lab*_p_hadron_lab+mass_p*mass_p) );

  /* 4-Vectors of colliding electron and proton in proton-at-rest frame */
  HLorentzVector *p4_electron_prest = new HLorentzVector( *p4_electron_lab );
  p4_electron_prest->boost(sqrt(p4_hadron_lab_invert->getQuare()),*p4_hadron_lab_invert);

  HLorentzVector *p4_hadron_prest = new HLorentzVector( *p4_hadron_lab );
  p4_hadron_prest->boost(sqrt(p4_hadron_lab_invert->getQuare()),*p4_hadron_lab_invert);

  if ( verbosity > -1 )
    {
      cout << "Electron and proton in laboratory frame:" << endl;
      p4_electron_lab->print();
      p4_hadron_lab->print();

      cout << "Electron and proton in proton-at-rest frame:" << endl;
      p4_electron_prest->print();
      p4_hadron_prest->print();
    }

  /* get instance of HepGenManager */
  _hgenManager = HGenManager::getInstance();
  _hgenManager->loadSettings( _datacardFile );
  _hgenManager->setupGenerator();

  /* set beam parameters */
  cout << "Colliding " << p4_electron_lab->getEnergy() << " GeV electron with " << p4_hadron_lab->getEnergy() << " GeV proton (laboratory frame)" << endl;
  cout << "----ELEPT (proton-at-rest): " << p4_electron_prest->getEnergy() << endl;
  _hgenManager->getParamManager()->getStruct()->ELEPT = p4_electron_prest->getEnergy();
  _hgenManager->getParamManager()->getStruct()->PARL.at(2) = p4_electron_prest->getEnergy();

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


int sHEPGen::process_event(PHCompositeNode *topNode) {

  if (verbosity > 1) cout << "sHEPGen::process_event - event: " << _eventcount << endl;

  _hgenManager->oneShot();

  if ( _detailed_debug )
    _hgenManager->getEvent()->printDebug();

  HEvent *evt_mc = _hgenManager->getEvent();

  /* Create HepMC GenEvent */
  HepMC::GenEvent* evt = new HepMC::GenEvent();

  /* define the units (Pythia uses GeV and mm) */
  evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

  /* add global information to the event */
  evt->set_event_number(_eventcount);

  /* Create single HepMC vertex for event - TODO: Allow for multiple vertices e.g. for decay particles */
  HepMC::GenVertex* hepmcvtx = new HepMC::GenVertex( HepMC::FourVector( 0,
                                                                        0,
                                                                        0,
                                                                        0 )
                                                     );

  /* Create HepMC particle records */
  HEventData* edata = evt_mc->getStruct();
  for ( unsigned p = 0; p < edata->listOfParticles.size(); p++ )
    {
      if (verbosity > 4)
        {
          cout << "______new particle_______" << endl;
          cout << "Index:  " << p+1 << endl;
          cout << "PID: " << edata->listOfParticles.at(p)->getParticleType() << " -- " << (edata->listOfParticles.at(p) == &(edata->incBeamParticle)) << endl;
          cout << "Particle aux flag: " << edata->listOfParticles.at(p)->getParticleAuxFlag() << endl;
          cout << "Particle origin: " << edata->listOfParticles.at(p)->getParticleOrigin() << endl;
          cout << "Particle daughter1: " << edata->listOfParticles.at(p)->getParticleDaughter1() << endl;
          cout << "Particle daughter2: " << edata->listOfParticles.at(p)->getParticleDaughter2() << endl;
        }

      HLorentzVector v4_particle_p = edata->listOfParticles.at(p)->getVector();
      HepMC::GenParticle *particle_hepmc = new HepMC::GenParticle( HepMC::FourVector(v4_particle_p.getVector().X(),
                                                                                     v4_particle_p.getVector().Y(),
                                                                                     v4_particle_p.getVector().Z(),
                                                                                     v4_particle_p.getEnergy()),
                                                                   edata->listOfParticles.at(p)->getParticleType() );
      particle_hepmc->set_status( edata->listOfParticles.at(p)->getParticleAuxFlag() );
      hepmcvtx->add_particle_out( particle_hepmc );
    }

  /* Add vertex to event */
  evt->add_vertex( hepmcvtx );

  /* pass HepMC to PHNode */
  bool success = _phhepmcevt->addEvent(evt);
  if (!success) {
    cout << "sHEPGen::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /* print outs */
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
