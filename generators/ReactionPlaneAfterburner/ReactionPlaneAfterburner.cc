
#include "ReactionPlaneAfterburner.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
//phhepmc
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

//hepmc
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>   // for GenParticle
#include <HepMC/GenVertex.h>     // for GenVertex, GenVertex::part...
#include <HepMC/SimpleVector.h>
#include <HepMC/HeavyIon.h>    // for HeavyIon
#pragma GCC diagnostic pop

#include <cassert>

//____________________________________________________________________________..
ReactionPlaneAfterburner::ReactionPlaneAfterburner(const std::string &name):
 SubsysReco(name)
{
   RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
}

//____________________________________________________________________________..
ReactionPlaneAfterburner::~ReactionPlaneAfterburner()
{
  gsl_rng_free(RandomGenerator);
}

//____________________________________________________________________________..
int ReactionPlaneAfterburner::Init(PHCompositeNode * /*topNode*/)
{
  unsigned int seed = PHRandomSeed();
  gsl_rng_set(RandomGenerator, seed);
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ReactionPlaneAfterburner::process_event(PHCompositeNode *topNode)
{
  //get event
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  for (auto & iter : *genevtmap)
  {
    PHHepMCGenEvent *genevt = iter.second;
    HepMC::GenEvent *evt = genevt->getEvent();
    assert(evt);
    HepMC::HeavyIon *hi = evt->heavy_ion();
    if (hi)
    {
      double psi = hi->event_plane_angle();
      if(psi != 0){
         std::cout << "ReactionPlaneAfterburner::process_event(PHCompositeNode *topNode) psi = " << psi << std::endl;
         std::cout<<"non-zero psi found, skipping"<<std::endl;
         return Fun4AllReturnCodes::EVENT_OK;
      }
      psi =  gsl_rng_uniform_pos(RandomGenerator) * 2 * M_PI;
      hi->set_event_plane_angle(psi);
      //rotate all particles and vertices
      for (HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p)
      {
        HepMC::FourVector v = (*p)->momentum();
        double px = v.px();
        double py = v.py();
        v.setPx(px * cos(psi) - py * sin(psi));
        v.setPy(px * sin(psi) + py * cos(psi));
        (*p)->set_momentum(v);
      }
      for (HepMC::GenEvent::vertex_iterator v = evt->vertices_begin(); v != evt->vertices_end(); ++v)
      {
        HepMC::FourVector pos = (*v)->position();
        double x = pos.x();
        double y = pos.y();
        pos.setX(x * cos(psi) - y * sin(psi));
        pos.setY(x * sin(psi) + y * cos(psi));
        (*v)->set_position(pos);
      }
    }
    else{
      std::cout<<"ReactionPlaneAfterburner::process_event: no heavy ion info found, exiting"<<std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }  
  return Fun4AllReturnCodes::EVENT_OK;
}




