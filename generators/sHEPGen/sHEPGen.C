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
  /* --> beam particle */
  HLorentzVector v4_beam = evt_mc->getBeam().getVector();
  int type_beam = evt_mc->getBeam().getParticleType();
  HepMC::GenParticle *particle_beam = new HepMC::GenParticle( HepMC::FourVector(v4_beam.getVector().X(),
										v4_beam.getVector().Y(),
										v4_beam.getVector().Z(),
										v4_beam.getEnergy()),
							      type_beam );
  particle_beam->set_status( 4 ); // status 4 = beam particle
  hepmcvtx->add_particle_in( particle_beam );

  /* --> target particle */
  int type_target = evt_mc->getRecoil().getParticleType();
  HepMC::GenParticle *particle_target = new HepMC::GenParticle( HepMC::FourVector(0,
										  0,
										  0,
										  evt_mc->getRecoil().getMass()),
								type_target );
  particle_target->set_status( 4 ); // status 4 = beam particle
  hepmcvtx->add_particle_in( particle_target );

  /* --> scattered beam particle */
  HLorentzVector v4_scat = evt_mc->getScat().getVector();
  int type_scat = evt_mc->getScat().getParticleType();
  HepMC::GenParticle *particle_scat = new HepMC::GenParticle( HepMC::FourVector(v4_scat.getVector().X(),
										v4_scat.getVector().Y(),
										v4_scat.getVector().Z(),
										v4_scat.getEnergy()),
							      type_scat );
  particle_scat->set_status( evt_mc->getScat().getParticleAuxFlag() );
  hepmcvtx->add_particle_out( particle_scat );

  /* --> recoil target particle */
  HLorentzVector v4_recoil = evt_mc->getRecoil().getVector();
  int type_recoil = evt_mc->getRecoil().getParticleType();
  HepMC::GenParticle *particle_recoil = new HepMC::GenParticle( HepMC::FourVector(v4_recoil.getVector().X(),
										v4_recoil.getVector().Y(),
										v4_recoil.getVector().Z(),
										v4_recoil.getEnergy()),
							      type_recoil );
  particle_recoil->set_status( evt_mc->getRecoil().getParticleAuxFlag() );
  hepmcvtx->add_particle_out( particle_recoil );

  /* --> output particle 1 */
  HLorentzVector v4_OutPart1 = evt_mc->getOutPart1().getVector();
  int type_OutPart1 = evt_mc->getOutPart1().getParticleType();
  if ( type_OutPart1 != 0 )
    {
      HepMC::GenParticle *particle_OutPart1 = new HepMC::GenParticle( HepMC::FourVector(v4_OutPart1.getVector().X(),
											v4_OutPart1.getVector().Y(),
											v4_OutPart1.getVector().Z(),
											v4_OutPart1.getEnergy()),
								      type_OutPart1 );
      particle_OutPart1->set_status( evt_mc->getOutPart1().getParticleAuxFlag() );
      hepmcvtx->add_particle_out( particle_OutPart1 );
    }

  /* --> output particle 2 */
  HLorentzVector v4_OutPart2 = evt_mc->getOutPart2().getVector();
  int type_OutPart2 = evt_mc->getOutPart2().getParticleType();
  if ( type_OutPart2 != 0 )
    {
      HepMC::GenParticle *particle_OutPart2 = new HepMC::GenParticle( HepMC::FourVector(v4_OutPart2.getVector().X(),
											v4_OutPart2.getVector().Y(),
											v4_OutPart2.getVector().Z(),
											v4_OutPart2.getEnergy()),
								      type_OutPart2 );
      particle_OutPart2->set_status( evt_mc->getOutPart2().getParticleAuxFlag() );
      hepmcvtx->add_particle_out( particle_OutPart2 );
    }

  /* --> output particle 3 */
  HLorentzVector v4_OutPart3 = evt_mc->getOutPart3().getVector();
  int type_OutPart3 = evt_mc->getOutPart3().getParticleType();
  if ( type_OutPart3 != 0 )
    {
      HepMC::GenParticle *particle_OutPart3 = new HepMC::GenParticle( HepMC::FourVector(v4_OutPart3.getVector().X(),
											v4_OutPart3.getVector().Y(),
											v4_OutPart3.getVector().Z(),
											v4_OutPart3.getEnergy()),
								      type_OutPart3 );
      particle_OutPart3->set_status( evt_mc->getOutPart3().getParticleAuxFlag() );
      hepmcvtx->add_particle_out( particle_OutPart3 );
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
