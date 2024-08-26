#include "TpcSeedFilter.h"

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/getClass.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

int TpcSeedFilter::process_event(PHCompositeNode* topNode)
{
  if (topNode == nullptr)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (Verbosity() > 1000)
  {
    topNode->print(); 
  }

  // Get the map of all the input seeds
  SvtxTrackMap* seeds = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!seeds)
  {
    std::cout << "Could not locate SvtxTrackMap node when running "
              << "\"TpcSeedFilter\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // iterate through the seeds and remove the ones that fail the cuts
  for (SvtxTrackMap::Iter track_pair = seeds->begin();
      track_pair != seeds->end(); ++track_pair)
  {
    auto& track = track_pair->second;
    if (track == nullptr) continue;
    
    if (
        ( _cut_on_min_pt && (track->get_pt() < _min_pt))
     || ( _cut_on_max_pt && (track->get_pt() > _max_pt))
     || ( _cut_on_min_eta && (track->get_eta() < _min_eta))
     || ( _cut_on_max_eta && (track->get_eta() > _max_eta))
     || ( track->size_cluster_keys() < _nclus_min )
     ) {
      if (Verbosity() > 4) {
        std::cout << " Cutting track: id(" << static_cast<int>(track_pair->first)
          << ")  nclusters: " // << std::static_cast<int>(track->size_cluster_keys() 
          << "   pt: " << track->get_pt() << "  eta: " << track->get_eta() << std::endl;
      }
      seeds->erase(track_pair->first);
      continue;
    }
    if (_must_span_sectors) {
      // check that there are clusters in at least 2 of the three layers of sectors
      bool in_0 = false;
      bool in_1 = false;
      bool in_2 = false;
      unsigned int sec_cnt = 0;
      for (auto key = track->begin_cluster_keys();
          key != track->end_cluster_keys();
          ++key)
      {
        unsigned int layer = TrkrDefs::getLayer(*key);
        if (Verbosity() > 4)
        {
          std::cout << ((int) layer) << " ";
        }

        if (layer < 23 && !in_0)
        {
          in_0 = true; 
          sec_cnt += 1;
        }
        else if (layer < 40 && !in_1)
        { 
          in_1 = true; 
          sec_cnt += 1;
        }
        else if ( layer < 49 && !in_2) // CHECK
        { 
          in_2 = true; 
          sec_cnt += 1;
        }

        if (sec_cnt >= _min_radial_sectors) break;
      }
      if (sec_cnt < _min_radial_sectors) {
        //cut the track
        seeds->erase(track_pair->first);

        if (Verbosity() > 4) {
          std::cout << " Cutting track: id(" << static_cast<int>(track_pair->first)
            << ") :  clusters only in ";
          if (in_0) std::cout << " inner ";
          if (in_1) std::cout << " middle ";
          if (in_2) std::cout << " outer ";
          std::cout << " radial sectors but needs at least " 
            << static_cast<int>(_min_radial_sectors) << std::endl;
        }
        continue;
      }
    }
  }
  if (Verbosity()>5) {
    std::cout << " Done cutting on QA for tracks." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
