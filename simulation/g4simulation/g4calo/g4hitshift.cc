#include "g4hitshift.h"

#include <g4detectors/PHG4CylinderGeom_Spacalv3.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                     // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <map>                                      // for operator!=, _Rb_t...
#include <string>
#include <utility>                                  // for pair

//____________________________________________________________________________..
g4hitshift::g4hitshift(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int g4hitshift::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_CEMC");
  if (hits)
  {
    //-----------------------------------------------------------------------
    // Loop over G4Hits
    //-----------------------------------------------------------------------
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      int scint_id = hit_iter->second->get_scint_id();

      PHG4CylinderGeom_Spacalv3::scint_id_coder decoder(scint_id);

      int sectornumber = decoder.sector_ID;
      int shiftedsector = sectornumber + 8;
      if (shiftedsector > 31)
      {
        shiftedsector = shiftedsector - 31;
      }
      PHG4CylinderGeom_Spacalv3::scint_id_coder encoder(shiftedsector, decoder.tower_ID, decoder.fiber_ID);

      hit_iter->second->set_scint_id(encoder.scint_ID);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
