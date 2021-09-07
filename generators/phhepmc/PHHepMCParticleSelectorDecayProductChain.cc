#include "PHHepMCParticleSelectorDecayProductChain.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>          // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>                 // for PHWHERE

#include <HepMC/GenEvent.h>              // for GenEvent::particle_const_ite...
#include <HepMC/GenParticle.h>           // for GenParticle
#include <HepMC/GenVertex.h>             // for GenVertex, GenVertex::partic...
#include <HepMC/IteratorRange.h>         // for ancestors, children, descend...
#include <HepMC/SimpleVector.h>          // for FourVector

#include <cstdlib>                      // for abs
#include <iostream>                      // for operator<<, basic_ostream::o...

class PHCompositeNode;

using namespace std;

PHHepMCParticleSelectorDecayProductChain::PHHepMCParticleSelectorDecayProductChain(const string& name)
  : SubsysReco(name)
  , _embedding_id(0)
{
  _theParticle = 11;
  return;
}

int PHHepMCParticleSelectorDecayProductChain::InitRun(PHCompositeNode* /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

HepMC::GenParticle* PHHepMCParticleSelectorDecayProductChain::GetParent(HepMC::GenParticle* p, HepMC::GenEvent* /*event*/)
{
  HepMC::GenParticle* parent = nullptr;
  if (!p->production_vertex()) return parent;

  for (HepMC::GenVertex::particle_iterator mother = p->production_vertex()->particles_begin(HepMC::ancestors);
       mother != p->production_vertex()->particles_end(HepMC::ancestors);
       ++mother)
  {
    for (unsigned int i = 0; i < _theAncestors.size(); i++)
    {
      if (abs((*mother)->pdg_id()) == _theAncestors[i])
      {
        parent = *mother;
        break;
      }
    }
    if (parent != nullptr) break;
  }

  return parent;
}

int PHHepMCParticleSelectorDecayProductChain::process_event(PHCompositeNode* topNode)
{
  if (_theParticle == 0 && _theDaughters.size() == 0)
  {
    cout << PHWHERE << "Doing nothing." << endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHHepMCGenEventMap* geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    cerr << "ERROR: PHHepMCGenEventMap node not found!" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHHepMCGenEvent* inEvent = geneventmap->get(_embedding_id);
  if (!inEvent)
  {
    cerr << "PHHepMCParticleSelectorDecayProductChain::process_event - WARNING: PHHepMCGenEvent with embedding ID " << _embedding_id << " is not found! Move on." << endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  HepMC::GenEvent* event = inEvent->getEvent();
  int npart = event->particles_size();
  int nvert = event->vertices_size();
  if (Verbosity() > 0) cout << "=========== Event " << event->event_number() << " contains " << npart << " particles and " << nvert << " vertices." << endl;

  // list of vertices to keep
  vector<HepMC::GenVertex> vkeep;

  if (_theParticle != 0)  // keep _theParticle and all its daughters (if any)
  {
    for (HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p)
    {
      // find _thePartcle
      if (abs((*p)->pdg_id()) == _theParticle)
      {
        // do we need to check for ancestors?
        if (_theAncestors.size() > 0)
        {
          HepMC::GenParticle* parent = GetParent(*p, event);
          if (parent)
          {
            vkeep.push_back(*(*p)->production_vertex());
          }
        }
        else
        {
          vkeep.push_back(*(*p)->production_vertex());
        }

        // do we need to keep the daughters?
        if (_theDaughters.size() > 0)
        {
          if ((*p)->end_vertex())
          {
            for (HepMC::GenVertex::particle_iterator des = (*p)->end_vertex()->particles_begin(HepMC::descendants);
                 des != (*p)->end_vertex()->particles_end(HepMC::descendants);
                 ++des)
            {
              for (unsigned int i = 0; i < _theDaughters.size(); i++)
              {
                if (abs((*des)->pdg_id()) == _theDaughters[i])
                {
                  vkeep.push_back(*(*p)->end_vertex());
                  break;
                }
              }
            }
          }
        }  // there are daughters

      }  // this is _theParticle
    }    // end loop over particles
  }
  else  // save only particles in _theDaughters list no matter where they came from
  {
    for (HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p)
    {
      for (unsigned int ip = 0; ip < _theDaughters.size(); ip++)
      {
        if (abs((*p)->pdg_id()) == _theDaughters[ip])
        {
          vkeep.push_back(*(*p)->production_vertex());
        }
      }
    }
  }

  // loop over vertices and keep only selected ones.
  for (HepMC::GenEvent::vertex_const_iterator v = event->vertices_begin(); v != event->vertices_end(); ++v)
  {
    bool goodvertex = false;
    for (unsigned int i = 0; i < vkeep.size(); i++)
    {
      HepMC::GenVertex tmp1 = (*(*v));
      HepMC::GenVertex tmp2 = vkeep[i];
      if (tmp1 == tmp2)
      {
        goodvertex = true;
      }
    }
    if (!goodvertex)
    {
      bool tmp = event->remove_vertex((*v));
      if (Verbosity() > 10 && tmp)
      {
        cout << PHWHERE << " Erasing empty vertex." << endl;
      }
    }
  }

  // clean up the vertices
  for (HepMC::GenEvent::vertex_const_iterator v = event->vertices_begin(); v != event->vertices_end(); ++v)
  {
    std::vector<HepMC::GenParticle*> removep;

    for (HepMC::GenVertex::particle_iterator itpart = (*v)->particles_begin(HepMC::children);
         itpart != (*v)->particles_end(HepMC::children);
         ++itpart)
    {
      bool keepparticle = false;
      if (abs((*itpart)->pdg_id()) == _theParticle)
      {
        keepparticle = true;
      }
      for (unsigned int j = 0; j < _theDaughters.size(); j++)
      {
        if (abs((*itpart)->pdg_id()) == _theDaughters[j] && (*itpart)->status() == 1)
        {
          keepparticle = true;
        }
      }
      if (!keepparticle)
      {
        removep.push_back((*itpart));
      }
    }  // end loop over particles in this vertex

    for (HepMC::GenVertex::particle_iterator itpart = (*v)->particles_begin(HepMC::parents);
         itpart != (*v)->particles_end(HepMC::parents);
         ++itpart)
    {
      bool keepparticle = false;
      if (abs((*itpart)->pdg_id()) == _theParticle)
      {
        keepparticle = true;
      }
      for (unsigned int j = 0; j < _theDaughters.size(); j++)
      {
        if (abs((*itpart)->pdg_id()) == _theDaughters[j] && (*itpart)->status() == 1)
        {
          keepparticle = true;
        }
      }
      if (!keepparticle)
      {
        removep.push_back((*itpart));
      }
    }  // end loop over particles in this vertex

    for (unsigned int k = 0; k < removep.size(); k++)
    {
      HepMC::GenParticle* tmp = (*v)->remove_particle(removep[k]);
      if (tmp->end_vertex())
      {
        delete tmp->end_vertex();
      }
    }
  }

  int partcount = 0;
  if (Verbosity() > 0)
  {
    cout << "FINAL Event " << event->event_number() << " contains " << event->particles_size() << " particles and " << event->vertices_size() << " vertices." << endl;
    cout << "FINAL LIST OF PARTICLES:" << endl;
  }
  for (HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p)
  {
    int pid = (*p)->pdg_id();
    int status = (*p)->status();
    double pz = ((*p)->momentum()).pz();
    double pt = ((*p)->momentum()).perp();
    double eta = ((*p)->momentum()).eta();
    double mass = ((*p)->momentum()).m();
    if (Verbosity() > 0)
    {
      cout << pid << " " << mass << " " << status << " " << pt << " " << pz << " " << eta << " " << (*p)->production_vertex() << " " << (*p)->end_vertex() << endl;
    }
    partcount++;
  }

  // if there is nothing to write out the code crashes
  if (partcount == 0)
  {
    if (Verbosity() > 0)
    {
      cout << "EVENT ABORTED: No particles to write out." << endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  else
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
}

void PHHepMCParticleSelectorDecayProductChain::SetParticle(const int pid)
{
  _theParticle = pid;
  return;
}

void PHHepMCParticleSelectorDecayProductChain::AddAncestor(const int pid)
{
  _theAncestors.push_back(pid);
  return;
}

void PHHepMCParticleSelectorDecayProductChain::AddDaughter(const int pid)
{
  _theDaughters.push_back(pid);
  return;
}
