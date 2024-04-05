#include "PHTrackCleaner.h"

#include "PHTrackCleaner.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>            // for cluskey, getLayer, TrkrId
#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>     // for sqrt, fabs, atan2, cos
#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for map
#include <set>       // for _Rb_tree_const_iterator
#include <utility>   // for pair, make_pair

//____________________________________________________________________________..
PHTrackCleaner::PHTrackCleaner(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
PHTrackCleaner::~PHTrackCleaner() = default;

//____________________________________________________________________________..
int PHTrackCleaner::InitRun(PHCompositeNode *topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  return ret;
}

//____________________________________________________________________________..
int PHTrackCleaner::process_event(PHCompositeNode * /*topNode*/)
{
  if (Verbosity() > 0)
  {
    std::cout << PHWHERE << " track map size " << _track_map->size() << std::endl;
  }

  std::set<unsigned int> track_keep_list;
  std::set<unsigned int> track_delete_list;

  unsigned int good_track = 0;  // for diagnostic output only
  unsigned int ok_track = 0;    // tracks to keep

  std::multimap<unsigned int, unsigned int> tpcid_track_mmap;
  std::set<unsigned int> tpc_id_set;
  // loop over the fitted tracks
  for (auto &it : *_track_map)
  {
    auto track_id = it.first;
    auto track = it.second;
    if (!track)
    {
      continue;
    }

    auto tpc_seed = track->get_tpc_seed();
    unsigned int tpc_index = _tpc_seed_map->find(tpc_seed);

    auto tpc_track_pair = std::make_pair(tpc_index, track_id);

    tpc_id_set.insert(tpc_index);
    tpcid_track_mmap.insert(tpc_track_pair);
  }

  if (Verbosity() > 0)
  {
    std::cout << " tpcid_track_mmap  size " << tpcid_track_mmap.size() << std::endl;
  }

  // loop over the TPC seed ID's

  for (unsigned int tpc_id : tpc_id_set)
  {
    if (Verbosity() > 1)
    {
      std::cout << " TPC ID " << tpc_id << std::endl;
    }

    auto tpc_range = tpcid_track_mmap.equal_range(tpc_id);

    unsigned int best_id = 99999;
    double min_chisq_df = 99999.0;
    unsigned int best_ndf = 1;
    for (auto it = tpc_range.first; it != tpc_range.second; ++it)
    {
      unsigned int track_id = it->second;

      // note that the track no longer exists if it failed in the Acts fitter
      _track = _track_map->get(track_id);

      if (_track)
      {
        if (_pp_mode)
        {
          // skip tracks with no assigned crossing number in pp mode
          if (_track->get_crossing() == SHRT_MAX)
          {
            if (Verbosity() > 0)
            {
              std::cout << "     skip  track ID " << track_id << " crossing " << _track->get_crossing() << " chisq " << _track->get_chisq()
                        << " ndf " << _track->get_ndf() << std::endl;
            }

            continue;
          }
        }

        // Find the remaining track with the best chisq/ndf

        unsigned int si_index = UINT_MAX;
        auto si_seed = _track->get_silicon_seed();
        if (si_seed)
        {
          si_index = _silicon_seed_map->find(si_seed);
        }

        if (Verbosity() > 1)
        {
          std::cout << "        track ID " << track_id << " tpc index " << tpc_id << " si index " << si_index << " crossing " << _track->get_crossing()
                    << " chisq " << _track->get_chisq() << " ndf " << _track->get_ndf() << " min_chisq_df " << min_chisq_df << std::endl;
        }

        // only accept tracks with ndf > min_ndf - very small ndf means something went wrong, as does ndf undefined
        if (_track->get_chisq() / _track->get_ndf() < min_chisq_df && _track->get_ndf() > min_ndf && _track->get_ndf() != UINT_MAX)
        {
          min_chisq_df = _track->get_chisq() / _track->get_ndf();
          best_id = track_id;
          best_ndf = _track->get_ndf();
        }
      }
    }

    if (best_id != 99999)
    {
      double qual = min_chisq_df;

      if (Verbosity() > 1)
      {
        std::cout << "        best track for tpc_id " << tpc_id << " has track_id " << best_id << " best_ndf " << best_ndf << " chisq/ndf " << qual << std::endl;
      }

      if (qual < quality_cut * 2)
      {
        track_keep_list.insert(best_id);
        ok_track++;
        if (qual < quality_cut)
        {
          good_track++;
        }
      }
    }
    else
    {
      if (Verbosity() > 1)
      {
        std::cout << "        no track exists  for tpc_id " << tpc_id << std::endl;
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << " Number of good tracks with qual < " << quality_cut << "  is " << good_track << " OK tracks " << ok_track << std::endl;
  }

  // make a list of tracks that did not make the keep list
  for (auto &track_it : *_track_map)
  {
    auto id = track_it.first;

    auto set_it = track_keep_list.find(id);
    if (set_it == track_keep_list.end())
    {
      if (Verbosity() > 1)
      {
        std::cout << "    add id " << id << " to track_delete_list " << std::endl;
      }
      track_delete_list.insert(id);
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << " track_delete_list size " << track_delete_list.size() << std::endl;
  }

  // delete failed tracks
  for (unsigned int it : track_delete_list)
  {
    if (Verbosity() > 1)
    {
      std::cout << " erasing track ID " << it << std::endl;
    }
    _track_map->erase(it);
  }

  if (Verbosity() > 0)
  {
    std::cout << "Track map size after choosing best silicon match: " << _track_map->size() << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackCleaner::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTrackCleaner::GetNodes(PHCompositeNode *topNode)
{
  _tpc_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (!_tpc_seed_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find TpcTrackSeedContainer: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _silicon_seed_map = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (!_silicon_seed_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SiliconTrackSeedContainer " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
