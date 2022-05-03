#include "PHSartre.h"
#include "PHSartreGenTrigger.h"

#include <sartre/Enumerations.h>            // for incoherent
#include <sartre/Event.h>                   // for Particle, Event
#include <sartre/EventGeneratorSettings.h>  // for EventGeneratorSettings


#include <phhepmc/PHHepMCGenHelper.h>       // for PHHepMCGenHelper

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>             // for SubsysReco

#include <sartre/Sartre.h>

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>                 // for TLorentzVector
#include <TParticlePDG.h>                   // for TParticlePDG

#include <CLHEP/Vector/LorentzVector.h>

#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>              // for GenParticle
#include <HepMC/GenVertex.h>                // for GenVertex
#include <HepMC/PdfInfo.h>                  // for PdfInfo
#include <HepMC/SimpleVector.h>             // for FourVector
#include <HepMC/Units.h>                    // for GEV, MM

#include <gsl/gsl_rng.h>                    // for gsl_rng_uniform

#include <cfloat>                           // for FLT_EPSILON
#include <cmath>                            // for M_PI
#include <cstdlib>                          // for getenv
#include <iostream>                         // for operator<<, endl, basic_o...
#include <memory>                           // for allocator_traits<>::value...

class PHHepMCGenEvent;

using namespace std;

PHSartre::PHSartre(const std::string &name)
  : SubsysReco(name)
{
  char *charPath = getenv("SARTRE_DIR");
  if (!charPath)
  {
    cout << "PHSartre::Could not find $SARTRE_DIR path!" << endl;
    return;
  }

  //
  //  Create the generator and initialize it.
  //  Once initialized you cannot (should not) change
  //  the settings w/o re-initialing sartre.
  //
  _sartre = new Sartre();

  PHHepMCGenHelper::set_embedding_id(1);  // default embedding ID to 1
}

PHSartre::~PHSartre()
{
  delete _sartre;
}

int PHSartre::Init(PHCompositeNode *topNode)
{
  if (!_configFile.empty())
  {
    bool ok = _sartre->init(_configFile);
    if (!ok)
    {
      cerr << "Initialization of sartre failed." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  else
  {
    cerr << "Sartre configuration file must be specified" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  settings = _sartre->runSettings();
  settings->list();

  create_node_tree(topNode);

  // event numbering will start from 1
  _eventcount = 0;

  decay = new TGenPhaseSpace();  // for VM decays
  daughterID = settings->userInt();
  if (daughterID && (settings->vectorMesonId() != 22))
  {
    doPerformDecay = true;
    daughterMasses[0] = settings->lookupPDG(daughterID)->Mass();
    daughterMasses[1] = settings->lookupPDG(-daughterID)->Mass();
    cout << "PHSartre: "
         << "Will decay vector meson: ";
    cout << "PHSartre: " << settings->lookupPDG(settings->vectorMesonId())->GetName();
    cout << "PHSartre: "
         << " -> ";
    cout << "PHSartre: " << settings->lookupPDG(daughterID)->GetName();
    cout << " ";
    cout << "PHSartre: " << settings->lookupPDG(-daughterID)->GetName();
    cout << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHSartre::End(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 1) cout << "PHSartre::End - I'm here!" << endl;

  cout << "PHSartre: "
       << " Total cross-section: " << _sartre->totalCrossSection() << " nb" << endl;
  _sartre->listStatus();

  cout << " *-------  Begin PHSARTRE Trigger Statistics  ----------------------"
       << "-------------------------------------------------* " << endl;
  cout << " |                                                                "
       << "                                                 | " << endl;
  cout << "                         PHSartre::End - " << _eventcount
       << " events passed trigger" << endl;
  cout << "                         Fraction passed: " << _eventcount
       << "/" << _gencount
       << " = " << _eventcount / float(_gencount) << endl;
  cout << " *-------  End PHSARTRE Trigger Statistics  ------------------------"
       << "-------------------------------------------------* " << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//-* print pythia config info
void PHSartre::print_config() const
{
  settings->list();
}

int PHSartre::process_event(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 1) cout << "PHSartre::process_event - event: " << _eventcount << endl;

  bool passedTrigger = false;
  Event *event = nullptr;

  TLorentzVector *eIn = nullptr;
  TLorentzVector *pIn = nullptr;
  TLorentzVector *eOut = nullptr;
  TLorentzVector *gamma = nullptr;
  TLorentzVector *vm = nullptr;
  TLorentzVector *PomOut = nullptr;
  TLorentzVector *pOut = nullptr;
  TLorentzVector *vmDecay1 = nullptr;
  TLorentzVector *vmDecay2 = nullptr;
  unsigned int preVMDecaySize = 0;

  while (!passedTrigger)
  {
    ++_gencount;

    // Generate a Sartre event
    event = _sartre->generateEvent();

    //
    //  If Sartre is run in UPC mode, half of the events needs to be
    //  rotated around and axis perpendicular to z:
    //  (only for symmetric events)
    //
    if (settings->UPC() and settings->A() == settings->UPCA())
    {
      randomlyReverseBeams(event);
    }

    // for sPHENIX/RHIC p+Au
    // (see comments in ReverseBeams)
    // reverse when the proton emits the virtual photon

    if (settings->UPC() and settings->A() == 197)
    {
      ReverseBeams(event);
    }

    // Set pointers to the parts of the event we will need:

    eIn = &event->particles[0].p;
    pIn = &event->particles[1].p;
    eOut = &event->particles[2].p;
    gamma = &event->particles[3].p;
    vm = &event->particles[4].p;
    PomOut = &event->particles[5].p;
    pOut = &event->particles[6].p;

    // To allow the triggering to work properly, we need to decay the vector meson here

    preVMDecaySize = event->particles.size();

    if (doPerformDecay)
    {
      if (decay->SetDecay(*vm, 2, daughterMasses))
      {
        double weight = decay->Generate();  // weight is always 1 here
        if ((weight - 1) > FLT_EPSILON)
        {
          cout << "PHSartre: Warning decay weight != 1, weight = " << weight << endl;
        }
        TLorentzVector *vmDaughter1 = decay->GetDecay(0);
        TLorentzVector *vmDaughter2 = decay->GetDecay(1);

        event->particles[4].status = 2;  // set VM status

        Particle vmDC1;
        vmDC1.index = event->particles.size();
        vmDC1.pdgId = daughterID;
        vmDC1.status = 1;  // final state
        vmDC1.p = *vmDaughter1;
        vmDC1.parents.push_back(4);
        event->particles.push_back(vmDC1);
        vmDecay1 = &event->particles[event->particles.size() - 1].p;

        Particle vmDC2;
        vmDC2.index = event->particles.size();
        vmDC2.pdgId = -daughterID;
        vmDC2.status = 1;  // final state
        vmDC2.p = *vmDaughter2;
        vmDC2.parents.push_back(4);
        event->particles.push_back(vmDC2);
        vmDecay2 = &event->particles[event->particles.size() - 1].p;
      }
      else
      {
        cout << "PHSartre: WARNING: Kinematics of Vector Meson does not allow decay!" << endl;
      }
    }

    // test trigger logic

    bool andScoreKeeper = true;
    if (Verbosity() > 2)
    {
      cout << "PHSartre::process_event - triggersize: " << _registeredTriggers.size() << endl;
    }

    for (unsigned int tr = 0; tr < _registeredTriggers.size(); tr++)
    {
      bool trigResult = _registeredTriggers[tr]->Apply(event);

      if (Verbosity() > 2)
      {
        cout << "PHSartre::process_event trigger: "
             << _registeredTriggers[tr]->GetName() << "  " << trigResult << endl;
      }

      if (_triggersOR && trigResult)
      {
        passedTrigger = true;
        break;
      }
      else if (_triggersAND)
      {
        andScoreKeeper &= trigResult;
      }

      if (Verbosity() > 2 && !passedTrigger)
      {
        cout << "PHSartre::process_event - failed trigger: "
             << _registeredTriggers[tr]->GetName() << endl;
      }
    }

    if ((andScoreKeeper && _triggersAND) || (_registeredTriggers.size() == 0))
    {
      passedTrigger = true;
    }
  }

  // fill HepMC object with event

  HepMC::GenEvent *genevent = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);

  // add some information to the event
  genevent->set_event_number(_eventcount);

  // Set the PDF information
  HepMC::PdfInfo pdfinfo;
  pdfinfo.set_scalePDF(event->Q2);
  genevent->set_pdf_info(pdfinfo);

  // We would also like to save:
  //
  // event->t;
  // event->x;
  // event->y;
  // event->s;
  // event->W;
  // event->xpom;
  // (event->polarization == transverse ? 0 : 1);
  // (event->diffractiveMode == coherent ? 0 : 1);
  //
  // but there doesn't seem to be a good place to do so
  // within the HepMC event information?
  //
  // t, W and Q^2 form a minial set of good variables for diffractive events
  // Maybe what I do is record the input particles to the event at the HepMC
  // vertices and reconstruct the kinematics from there?

  // Create HepMC vertices and add final state particles to them

  // First, the emitter(electron)-virtual photon vertex:

  HepMC::GenVertex *egammavtx = new HepMC::GenVertex(CLHEP::HepLorentzVector(0.0, 0.0, 0.0, 0.0));
  genevent->add_vertex(egammavtx);

  egammavtx->add_particle_in(
      new HepMC::GenParticle(CLHEP::HepLorentzVector(eIn->Px(),
                                                     eIn->Py(),
                                                     eIn->Pz(),
                                                     eIn->E()),
                             event->particles[0].pdgId,
                             3));

  HepMC::GenParticle *hgamma = new HepMC::GenParticle(CLHEP::HepLorentzVector(gamma->Px(),
                                                                              gamma->Py(),
                                                                              gamma->Pz(),
                                                                              gamma->E()),
                                                      event->particles[3].pdgId,
                                                      3);

  egammavtx->add_particle_out(hgamma);

  egammavtx->add_particle_out(
      new HepMC::GenParticle(CLHEP::HepLorentzVector(eOut->Px(),
                                                     eOut->Py(),
                                                     eOut->Pz(),
                                                     eOut->E()),
                             event->particles[2].pdgId,
                             1));

  // Next, the hadron-pomeron vertex:

  HepMC::GenVertex *ppomvtx = new HepMC::GenVertex(CLHEP::HepLorentzVector(0.0, 0.0, 0.0, 0.0));
  genevent->add_vertex(ppomvtx);

  ppomvtx->add_particle_in(
      new HepMC::GenParticle(CLHEP::HepLorentzVector(pIn->Px(),
                                                     pIn->Py(),
                                                     pIn->Pz(),
                                                     pIn->E()),
                             event->particles[1].pdgId,
                             3));

  HepMC::GenParticle *hPomOut = new HepMC::GenParticle(CLHEP::HepLorentzVector(PomOut->Px(),
                                                                               PomOut->Py(),
                                                                               PomOut->Pz(),
                                                                               PomOut->E()),
                                                       event->particles[5].pdgId,
                                                       3);

  ppomvtx->add_particle_out(hPomOut);

  // If this is a nuclear breakup, add in the nuclear fragments
  // Otherwise, add in the outgoing hadron

  //If the event is incoherent, and nuclear breakup is enabled, fill the remnants to the tree
  if (settings->enableNuclearBreakup() and event->diffractiveMode == incoherent)
  {
    for (unsigned int iParticle = 7; iParticle < preVMDecaySize; iParticle++)
    {
      if (event->particles[iParticle].status == 1)
      {  // Final-state particle
        const Particle &particle = event->particles[iParticle];
        ppomvtx->add_particle_out(
            new HepMC::GenParticle(CLHEP::HepLorentzVector(particle.p.Px(),
                                                           particle.p.Py(),
                                                           particle.p.Pz(),
                                                           particle.p.E()),
                                   particle.pdgId,
                                   1));
      }
    }
  }
  else
  {
    ppomvtx->add_particle_out(
        new HepMC::GenParticle(CLHEP::HepLorentzVector(pOut->Px(),
                                                       pOut->Py(),
                                                       pOut->Pz(),
                                                       pOut->E()),
                               event->particles[6].pdgId,
                               1));
  }

  // The Pomeron-Photon vertex

  HepMC::GenVertex *gammapomvtx = new HepMC::GenVertex(CLHEP::HepLorentzVector(0.0, 0.0, 0.0, 0.0));
  genevent->add_vertex(gammapomvtx);

  gammapomvtx->add_particle_in(hgamma);
  gammapomvtx->add_particle_in(hPomOut);

  int isVMFinal = 1;
  if (doPerformDecay) isVMFinal = 2;

  HepMC::GenParticle *hvm = new HepMC::GenParticle(CLHEP::HepLorentzVector(vm->Px(),
                                                                           vm->Py(),
                                                                           vm->Pz(),
                                                                           vm->E()),
                                                   event->particles[4].pdgId,
                                                   isVMFinal);

  gammapomvtx->add_particle_out(hvm);

  // Add the VM decay to the event

  if (doPerformDecay)
  {
    if (vmDecay1 && vmDecay2)
    {
      HepMC::GenVertex *fvtx = new HepMC::GenVertex(CLHEP::HepLorentzVector(0.0, 0.0, 0.0, 0.0));
      genevent->add_vertex(fvtx);

      fvtx->add_particle_in(hvm);

      fvtx->add_particle_out(
          new HepMC::GenParticle(CLHEP::HepLorentzVector(vmDecay1->Px(),
                                                         vmDecay1->Py(),
                                                         vmDecay1->Pz(),
                                                         vmDecay1->E()),
                                 daughterID,
                                 1));
      fvtx->add_particle_out(
          new HepMC::GenParticle(CLHEP::HepLorentzVector(vmDecay2->Px(),
                                                         vmDecay2->Py(),
                                                         vmDecay2->Pz(),
                                                         vmDecay2->E()),
                                 -daughterID,
                                 1));
    }
    else
    {
      cout << "PHSartre: WARNING: Kinematics of Vector Meson does not allow decay!" << endl;
    }
  }

  // pass HepMC to PHNode

  PHHepMCGenEvent *success = PHHepMCGenHelper::insert_event(genevent);
  if (!success)
  {
    cout << "PHSartre::process_event - Failed to add event to HepMC record!" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // print outs

  if (Verbosity() > 2) cout << "PHSartre::process_event - FINISHED WHOLE EVENT" << endl;

  ++_eventcount;
  return Fun4AllReturnCodes::EVENT_OK;
}



int PHSartre::ResetEvent(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHSartre::register_trigger(PHSartreGenTrigger *theTrigger)
{
  if (Verbosity() > 1) cout << "PHSartre::registerTrigger - trigger " << theTrigger->GetName() << " registered" << endl;
  _registeredTriggers.push_back(theTrigger);
}

// UPC only
void PHSartre::randomlyReverseBeams(Event *myEvent)
{
  if (gsl_rng_uniform(PHHepMCGenHelper::get_random_generator()) > 0.5)
  {
    for (unsigned int i = 0; i < myEvent->particles.size(); i++)
      myEvent->particles.at(i).p.RotateX(M_PI);
  }
}

// Used to rotate into the sPHENIX/RHIC frame
// The photon emitting beam is always pz<0, and in sPHENIX this will be
// the ion direction. So what you want to run are two files with:
// A=1 UPCA=197 (no reversal, Au ion emits photon)
// A=197 UPCA=1 (reversal required, proton emits photon)
void PHSartre::ReverseBeams(Event *myEvent)
{
  for (unsigned int i = 0; i < myEvent->particles.size(); i++)
    myEvent->particles.at(i).p.RotateX(M_PI);
}
