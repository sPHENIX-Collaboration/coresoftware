#include "g4hitshifthcal.h"

// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for hit_idbits

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <cassert>
#include <sstream>
#include <string>

//____________________________________________________________________________..
g4hitshifthcal::g4hitshifthcal(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int g4hitshifthcal::process_event(PHCompositeNode *topNode)
{
  PHG4HitContainer *hitsin = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALIN");
  if (hitsin)
  {
    //-----------------------------------------------------------------------
    // Loop over G4Hits
    //-----------------------------------------------------------------------
    PHG4HitContainer::ConstRange hit_range = hitsin->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      int introw = (hit_iter->second->get_hit_id() >> PHG4HitDefs::hit_idbits);

      // Get the original hit_id
      PHG4HitDefs::keytype original_hit_id = hit_iter->second->get_hit_id();

      // Get the lowBits from the original hit_id
      PHG4HitDefs::keytype lowBits = original_hit_id & ((1ULL << PHG4HitDefs::hit_idbits) - 1);
      ;

      // shift row up by 4
      int newrow = introw + 4;

      if (newrow >= 256)
      {
        newrow = newrow - 256;
      }
      // Clear the high-order bits of hit_id and set them to the new value
      PHG4HitDefs::keytype new_hit_id = (static_cast<PHG4HitDefs::keytype>(newrow) << PHG4HitDefs::hit_idbits) | lowBits;

      hit_iter->second->set_hit_id(new_hit_id);
    }
  }
  // for ohcal
  PHG4HitContainer *hitsout = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_HCALOUT");
  if (hitsout)
  {
    //-----------------------------------------------------------------------
    // Loop over G4Hits
    //-----------------------------------------------------------------------
    PHG4HitContainer::ConstRange hit_range = hitsout->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      int introw = (hit_iter->second->get_hit_id() >> PHG4HitDefs::hit_idbits);

      // Get the original hit_id
      PHG4HitDefs::keytype original_hit_id = hit_iter->second->get_hit_id();

      // Get the lowBits from the original hit_id
      PHG4HitDefs::keytype lowBits = original_hit_id & ((1ULL << PHG4HitDefs::hit_idbits) - 1);

      // shift row up by 5
      int newrow = introw + 5;

      if (newrow >= 320)
      {
        newrow = newrow - 320;
      }
      // Clear the high-order bits of hit_id and set them to the new value
      PHG4HitDefs::keytype new_hit_id = (static_cast<PHG4HitDefs::keytype>(newrow) << PHG4HitDefs::hit_idbits) | lowBits;

      hit_iter->second->set_hit_id(new_hit_id);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
