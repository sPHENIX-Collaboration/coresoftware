#include "sHEPGen.h"

#include <phhepmc/PHHepMCGenEvent.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>

/* HEPGen includes */
#include <hgenmanager.h>

/* HepMC includes */
#include <HepMC/GenEvent.h>

using namespace std;

typedef PHIODataNode<PHObject> PHObjectNode_t;

sHEPGen::sHEPGen(const std::string &name)
  : SubsysReco(name)
  , _eventcount(0)
  , _p_electron_lab(-20)
  , _p_hadron_lab(250)
  , _p4_electron_lab(nullptr)
  , _p4_hadron_lab(nullptr)
  , _p4_hadron_lab_invert(nullptr)
  , _p4_electron_prest(nullptr)
  , _p4_hadron_prest(nullptr)
  , _hgenManager(NULL)
  , _datacardFile("hepgen_dvcs.data")
{
  PHHepMCGenHelper::set_embedding_id(1);  // default embedding ID to 1
}

sHEPGen::~sHEPGen()
{
}

int sHEPGen::Init(PHCompositeNode *topNode)
{
  printlogo();

  /* electron and proton mass */
  double mass_e = 5.109989e-4;
  double mass_p = 9.382720e-1;

  /* 4-Vectors of colliding electron and proton in laboratory frame */
  _p4_electron_lab = new HLorentzVector(0.,
                                        0.,
                                        _p_electron_lab,
                                        sqrt(_p_electron_lab * _p_electron_lab + mass_e * mass_e));

  _p4_hadron_lab = new HLorentzVector(0.,
                                      0.,
                                      _p_hadron_lab,
                                      sqrt(_p_hadron_lab * _p_hadron_lab + mass_p * mass_p));

  /* The current version of HEPGen supports only a "Fixed Target" mode, i.e.
     the target (proton) is assumed at rest. Therefore, need to boost collision
     from the laboratory frame to the proton-at-rest frame. */
  _p4_hadron_lab_invert = new HLorentzVector(0.,
                                             0.,
                                             -1 * _p_hadron_lab,
                                             sqrt(_p_hadron_lab * _p_hadron_lab + mass_p * mass_p));

  /* 4-Vectors of colliding electron and proton in proton-at-rest frame */
  _p4_electron_prest = new HLorentzVector(*_p4_electron_lab);
  _p4_electron_prest->boost(sqrt(_p4_hadron_lab_invert->getQuare()), *_p4_hadron_lab_invert);

  _p4_hadron_prest = new HLorentzVector(*_p4_hadron_lab);
  _p4_hadron_prest->boost(sqrt(_p4_hadron_lab_invert->getQuare()), *_p4_hadron_lab_invert);

  if (verbosity > 2)
  {
    cout << "Electron and proton in laboratory frame:" << endl;
    _p4_electron_lab->print();
    _p4_hadron_lab->print();

    cout << "Electron and proton in proton-at-rest frame:" << endl;
    _p4_electron_prest->print();
    _p4_hadron_prest->print();
  }

  /* flip sign of electron momentum z-component for HEPGen */
  HVector3 flip_pz(_p4_electron_prest->getVector());
  flip_pz.setZ(flip_pz.Z() * -1);
  _p4_electron_prest->setVector(flip_pz);

  /* get instance of HepGenManager */
  _hgenManager = HGenManager::getInstance();
  _hgenManager->loadSettings(_datacardFile);
  _hgenManager->setupGenerator();

  /* set beam parameters */
  cout << "Colliding " << _p4_electron_lab->getVector().Z() << " GeV electron with " << _p4_hadron_lab->getVector().Z() << " GeV proton (laboratory frame)" << endl;
  cout << "----ELEPT (proton-at-rest): " << _p4_electron_prest->getVector().Z() << " GeV" << endl;
  _hgenManager->getParamManager()->getStruct()->ELEPT = _p4_electron_prest->getVector().Z();
  _hgenManager->getParamManager()->getStruct()->PARL.at(2) = _p4_electron_prest->getVector().Z();

  /* set random seed */
  unsigned int seed = PHRandomSeed();
  _hgenManager->setSeed(seed);

  /* enable detailed event record printput for debugging */
  if (verbosity > 2)
    _hgenManager->enableDebug();

  create_node_tree(topNode);

  _eventcount = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

int sHEPGen::End(PHCompositeNode *topNode)
{
  cout << "Reached the sHEPGen::End()" << endl;

  //-* dump out closing info (cross-sections, etc)

  if (verbosity > 1) cout << "sHEPGen::End - I'm here!" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int sHEPGen::process_event(PHCompositeNode *topNode)
{
  if (verbosity > 1) cout << "sHEPGen::process_event - event: " << _eventcount << endl;

  _hgenManager->oneShot();

  if (verbosity > 2)
    _hgenManager->getEvent()->printDebug();

  HEvent *evt_mc = _hgenManager->getEvent();

  /* Create HepMC GenEvent */
  HepMC::GenEvent *evt = new HepMC::GenEvent();

  /* define the units (Pythia uses GeV and mm) */
  evt->use_units(HepMC::Units::GEV, HepMC::Units::MM);

  /* add global information to the event */
  evt->set_event_number(_eventcount);

  /* Set the PDF information */
  HepMC::PdfInfo pdfinfo;
  pdfinfo.set_scalePDF(evt_mc->getQsq());
  pdfinfo.set_x2(evt_mc->getXbj());
  evt->set_pdf_info(pdfinfo);

  /* Set additional event parameters */
  evt->set_event_scale(evt_mc->getQsq());

  /* Event kinematics */
  /* @TODO How can we store this information in HepMC record? */
  //evt_mc->getNu();
  //evt_mc->getY();
  //evt_mc->getWsq();
  //evt_mc->getXbj();
  //evt_mc->getT();
  //evt_mc->getQsq();
  //evt_mc->getS();
  //evt_mc->getEmiss();

  /* Create single HepMC vertex for event */
  /* @TODO: Implement multiple vertices e.g. for decay particles */
  HepMC::GenVertex *hepmcvtx = new HepMC::GenVertex(HepMC::FourVector(0,
                                                                      0,
                                                                      0,
                                                                      0));

  /* Create HepMC particle records */
  HEventData *edata = evt_mc->getStruct();
  for (unsigned p = 0; p < edata->listOfParticles.size(); p++)
  {
    if (verbosity > 4)
    {
      cout << "______new particle_______" << endl;
      cout << "Index:  " << p + 1 << endl;
      cout << "PID: " << edata->listOfParticles.at(p)->getParticleType() << " -- " << (edata->listOfParticles.at(p) == &(edata->incBeamParticle)) << endl;
      cout << "Particle aux flag: " << edata->listOfParticles.at(p)->getParticleAuxFlag() << endl;
      cout << "Particle origin: " << edata->listOfParticles.at(p)->getParticleOrigin() << endl;
      cout << "Particle daughter1: " << edata->listOfParticles.at(p)->getParticleDaughter1() << endl;
      cout << "Particle daughter2: " << edata->listOfParticles.at(p)->getParticleDaughter2() << endl;
    }

    HLorentzVector v4_particle_p = edata->listOfParticles.at(p)->getVector();

    /* flip z axis component of particle momentum: sPHENIX EIC detector coordinate system uses
         electron flying in 'negative z' direction */
    HVector3 flip_pz(v4_particle_p.getVector());
    flip_pz.setZ(flip_pz.Z() * -1);
    v4_particle_p.setVector(flip_pz);

    /* Boost particle from proton-at-rest frame to laboratory frame */
    HLorentzVector v4_particle_p_lab(v4_particle_p);

    v4_particle_p_lab.boost(sqrt(_p4_hadron_lab->getQuare()), *_p4_hadron_lab);

    if (verbosity > 2)
    {
      cout << "EVENT RECORD particle: " << edata->listOfParticles.at(p)->getParticleType()
           << " (status: " << edata->listOfParticles.at(p)->getParticleAuxFlag() << ")" << endl;
      cout << " --> PROTON REST frame: ";
      v4_particle_p.print();
      cout << " --> LABORATORY  frame: ";
      v4_particle_p_lab.print();
    }

    HepMC::GenParticle *particle_hepmc = new HepMC::GenParticle(HepMC::FourVector(v4_particle_p_lab.getVector().X(),
                                                                                  v4_particle_p_lab.getVector().Y(),
                                                                                  v4_particle_p_lab.getVector().Z(),
                                                                                  v4_particle_p_lab.getEnergy()),
                                                                edata->listOfParticles.at(p)->getParticleType());
    particle_hepmc->set_status(edata->listOfParticles.at(p)->getParticleAuxFlag());
    hepmcvtx->add_particle_out(particle_hepmc);
  }

  /* Add vertex to event */
  evt->add_vertex(hepmcvtx);

  /* pass HepMC to PHNode */
  PHHepMCGenEvent *success = PHHepMCGenHelper::insert_event(evt);
  if (!success)
  {
    cout << "sHEPGen::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /* print outs */
  if (verbosity > 2) cout << "sHEPGen::process_event - FINISHED WHOLE EVENT" << endl;

  ++_eventcount;

  return Fun4AllReturnCodes::EVENT_OK;
}

void sHEPGen::printlogo()
{
  cout << endl
       << endl;
  for (int i = 0; i < hepconst::logolength; i++)
    cout << hepconst::logo[i] << endl;
  cout << endl
       << endl;
}
