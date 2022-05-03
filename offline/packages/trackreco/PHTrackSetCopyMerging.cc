#include "PHTrackSetCopyMerging.h"

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                      // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>                         // for PHWHERE

#include <iostream>                              // for operator<<, basic_os...

using namespace std;

PHTrackSetCopyMerging::PHTrackSetCopyMerging(const std::string& name)
  : PHTrackSetMerging(name)
{
}


int PHTrackSetCopyMerging::Process(PHCompositeNode* /*topNode*/)
{

  if (!_track_map_in1){
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name_in1  << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!_track_map_in2){
    cerr << PHWHERE << " ERROR: Can't find " << _track_map_name_in2  << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  //copy track map 1 to output track map
  /*
  // make a list of tracks that did not make the keep list
  for(auto track_it = _track_map_in1->begin(); track_it != _track_map_in1->end(); ++track_it)
    {
      // auto id = track_it->first;
      SvtxTrack *_track = track_it->second;

      _track_map_out->insert(_track);
  */
  //copy track map 2 to output track map
  for(auto track_it = _track_map_in2->begin(); track_it != _track_map_in2->end(); ++track_it)
    {
      //auto id = track_it->first;
      SvtxTrack *_track = track_it->second;

      //      _track_map_out->insert(_track);
      _track_map_in1->insert(_track);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
/*
int PHTrackSetCopyMerging::Setup(PHCompositeNode* topNode)
{
  int ret = CreateNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackSetMerging::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
*/
