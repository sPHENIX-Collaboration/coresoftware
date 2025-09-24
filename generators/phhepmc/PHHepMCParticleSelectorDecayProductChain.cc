#include "PHHepMCParticleSelectorDecayProductChain.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <HepMC/GenEvent.h>       // for GenEvent::particle_const_ite...
#include <HepMC/GenParticle.h>    // for GenParticle
#include <HepMC/GenVertex.h>      // for GenVertex, GenVertex::partic...
#include <HepMC/IteratorRange.h>  // for ancestors, children, descend...
#include <HepMC/SimpleVector.h>   // for FourVector

#include <cstdlib>   // for abs
#include <iostream>  // for operator<<, basic_ostream::o...

class PHCompositeNode;

PHHepMCParticleSelectorDecayProductChain::PHHepMCParticleSelectorDecayProductChain(const std::string& name)
  : SubsysReco(name)
{
  return;
}

HepMC::GenParticle* PHHepMCParticleSelectorDecayProductChain::GetParent(HepMC::GenParticle* p, HepMC::GenEvent* /*event*/)
{
  HepMC::GenParticle* parent = nullptr;
  if (!p->production_vertex())
  {
    return parent;
  }

  for (HepMC::GenVertex::particle_iterator mother = p->production_vertex()->particles_begin(HepMC::ancestors);
       mother != p->production_vertex()->particles_end(HepMC::ancestors);
       ++mother)
  {
    for (int _theAncestor : _theAncestors)
    {
      if (abs((*mother)->pdg_id()) == _theAncestor)
      {
        parent = *mother;
        break;
      }
    }
    if (parent != nullptr)
    {
      break;
    }
  }

  return parent;
}

int PHHepMCParticleSelectorDecayProductChain::process_event(PHCompositeNode* topNode)
{
  if (_theParticle == 0 && _theDaughters.empty())
  {
    std::cout << PHWHERE << "Doing nothing." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  PHHepMCGenEventMap* geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!geneventmap)
  {
    std::cout << "ERROR: PHHepMCGenEventMap node not found!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHHepMCGenEvent* inEvent = geneventmap->get(_embedding_id);
  if (!inEvent)
  {
    std::cout << "PHHepMCParticleSelectorDecayProductChain::process_event - WARNING: PHHepMCGenEvent with embedding ID " << _embedding_id << " is not found! Move on." << std::endl;
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  HepMC::GenEvent* event = inEvent->getEvent();
  int npart = event->particles_size();
  int nvert = event->vertices_size();
  if (Verbosity() > 0)
  {
    std::cout << "=========== Event " << event->event_number() << " contains " << npart << " particles and " << nvert << " vertices." << std::endl;
  }

  // list of vertices to keep
  std::vector<HepMC::GenVertex> vkeep;

  if (_theParticle != 0)  // keep _theParticle and all its daughters (if any)
  {
    for (HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p)
    {
      // find _thePartcle
      if (abs((*p)->pdg_id()) == _theParticle)
      {
        // do we need to check for ancestors?
        if (!_theAncestors.empty())
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
        if (!_theDaughters.empty())
        {
          if ((*p)->end_vertex())
          {
            for (HepMC::GenVertex::particle_iterator des = (*p)->end_vertex()->particles_begin(HepMC::descendants);
                 des != (*p)->end_vertex()->particles_end(HepMC::descendants);
                 ++des)
            {
              for (int _theDaughter : _theDaughters)
              {
                if (abs((*des)->pdg_id()) == _theDaughter)
                {
                  vkeep.push_back(*(*p)->end_vertex());
                  break;
                }
              }
            }
          }
        }  // there are daughters

      }  // this is _theParticle
    }  // end loop over particles
  }
  else  // save only particles in _theDaughters list no matter where they came from
  {
    for (HepMC::GenEvent::particle_const_iterator p = event->particles_begin(); p != event->particles_end(); ++p)
    {
      for (int _theDaughter : _theDaughters)
      {
        if (abs((*p)->pdg_id()) == _theDaughter)
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
    for (const auto& tmp2 : vkeep)
    {
      HepMC::GenVertex tmp1 = (*(*v));
      if (tmp1 == tmp2)
      {
        goodvertex = true;
        break;
      }
    }
    if (!goodvertex)
    {
      bool tmp = event->remove_vertex((*v));
      if (Verbosity() > 10 && tmp)
      {
        std::cout << PHWHERE << " Erasing empty vertex." << std::endl;
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
      for (int _theDaughter : _theDaughters)
      {
        if (abs((*itpart)->pdg_id()) == _theDaughter && (*itpart)->status() == 1)
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
      for (int _theDaughter : _theDaughters)
      {
        if (abs((*itpart)->pdg_id()) == _theDaughter && (*itpart)->status() == 1)
        {
          keepparticle = true;
        }
      }
      if (!keepparticle)
      {
        removep.push_back((*itpart));
      }
    }  // end loop over particles in this vertex

    for (auto& k : removep)
    {
      HepMC::GenParticle* tmp = (*v)->remove_particle(k);
      if (tmp->end_vertex())
      {
        delete tmp->end_vertex();
      }
    }
  }

  int partcount = 0;
  if (Verbosity() > 0)
  {
    std::cout << "FINAL Event " << event->event_number() << " contains " << event->particles_size() << " particles and " << event->vertices_size() << " vertices." << std::endl;
    std::cout << "FINAL LIST OF PARTICLES:" << std::endl;
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
      std::cout << pid << " " << mass << " " << status << " " << pt << " " << pz << " " << eta << " " << (*p)->production_vertex() << " " << (*p)->end_vertex() << std::endl;
    }
    partcount++;
  }

  // if there is nothing to write out the code crashes
  if (partcount == 0)
  {
    if (Verbosity() > 0)
    {
      std::cout << "EVENT ABORTED: No particles to write out." << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
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
