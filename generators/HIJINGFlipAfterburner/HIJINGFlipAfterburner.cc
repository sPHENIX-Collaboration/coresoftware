
#include "HIJINGFlipAfterburner.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
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
#pragma GCC diagnostic pop



#include <cassert>

//____________________________________________________________________________..
HIJINGFlipAfterburner::HIJINGFlipAfterburner(const std::string &name) : SubsysReco(name)
{
}

//____________________________________________________________________________..
HIJINGFlipAfterburner::~HIJINGFlipAfterburner()
{
}


//____________________________________________________________________________..
int HIJINGFlipAfterburner::process_event(PHCompositeNode *topNode)
{
  if (!doFlip)
  {
    doFlip = true;
  }
  else
  {
    doFlip = false;
    PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
    for (PHHepMCGenEventMap::Iter iter = genevtmap->begin(); iter != genevtmap->end(); ++iter)
    {
      PHHepMCGenEvent *genevt = iter->second;
      HepMC::GenEvent *evt = genevt->getEvent();
      assert(evt);
      flipZDirection(evt);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}



void HIJINGFlipAfterburner::flipZDirection(HepMC::GenEvent *event)
{
  assert(event);

  for (auto v = event->vertices_begin(); v != event->vertices_end(); ++v)
  {
    HepMC::FourVector position = (*v)->position();

    position.setZ(-position.z());
    position.setX(-position.x());

    (*v)->set_position(position);
  }

  for (auto p = event->particles_begin(); p != event->particles_end(); ++p)
  {
    HepMC::FourVector momentum = (*p)->momentum();

    momentum.setPz(-momentum.pz());
    momentum.setPx(-momentum.px());

    (*p)->set_momentum(momentum);
  }
}
