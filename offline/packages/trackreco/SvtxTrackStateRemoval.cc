
#include "SvtxTrackStateRemoval.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>

//____________________________________________________________________________..
SvtxTrackStateRemoval::SvtxTrackStateRemoval(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
SvtxTrackStateRemoval::~SvtxTrackStateRemoval()
{
}

//____________________________________________________________________________..
int SvtxTrackStateRemoval::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTrackStateRemoval::InitRun(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTrackStateRemoval::process_event(PHCompositeNode* topNode)
{
  auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
  {
    std::cout << PHWHERE << "No track map on node tree, can't continue."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto tpotgeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!tpotgeom)
  {
    std::cout << PHWHERE << "No micromegas geometry, can't continue."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /// get the last tracking layer
  auto layergeom = tpotgeom->GetLayerGeom(56);

  const float lastradius = layergeom->get_radius();
  const float lastthickness = layergeom->get_thickness();
  const float lasttrackingradius = lastradius + lastthickness / 2.;

  for (auto& [key, track] : *trackmap)
  {
    for (auto iter = track->begin_states(); iter != track->end_states(); ++iter)
    {
      /// Don't erase the PCA state information
      if (iter == track->begin_states())
      {
        continue;
      }

      float pathlength = iter->second->get_pathlength();
      if (pathlength < lasttrackingradius)
      {
        track->erase_state(pathlength);
      }
    }

    if (Verbosity() > 1)
    {
      track->identify();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int SvtxTrackStateRemoval::End(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
