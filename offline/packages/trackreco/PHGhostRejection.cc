#include "PHGhostRejection.h"

#include "PHGhostRejection.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for cluskey, getLayer, TrkrId

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <cmath>     // for sqrt, fabs, atan2, cos
#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for map
#include <set>       // for _Rb_tree_const_iterator
#include <utility>   // for pair, make_pair

//____________________________________________________________________________..
bool PHGhostRejection::cut_from_clusters(int itrack) {
  const auto& track = seeds[itrack];
  unsigned int nclusters = track.size_cluster_keys();
  if (nclusters < _min_clusters)
  {
    m_rejected[itrack] = true;
    if (m_verbosity > 3) {
      std::cout << " rejecting track ID (" << ((int)itrack)
        <<") because n-clusters(" << ((int)nclusters) <<")" << std::endl;
    }
    return true;
  }
  if (_must_span_sectors)
  {
    // check that there are clusters in at least 2 of the three layers of sectors
    bool in_two_sectors = false;
    bool in_0 = false;
    bool in_1 = false;
    bool in_2 = false;
    if (m_verbosity > 3)
    {
      std::cout << " layers in track: ";
    }
    for (auto key = track.begin_cluster_keys();
        key != track.end_cluster_keys();
        ++key)
    {
      unsigned int layer = TrkrDefs::getLayer(*key);
      if (m_verbosity > 4)
      {
        std::cout << ((int) layer) << " ";
      }

      if (layer < 23)
      { in_0 = true; }
      else if (layer < 40)
      { in_1 = true; }
      else
      { in_2 = true; }

      if ((in_0+in_1+in_2)>1)
      {
        in_two_sectors = true;
        break;
      }
    }
    if (m_verbosity > 3)
    { std::cout << std::endl; }
    if (!in_two_sectors)
    {
      m_rejected[itrack] = true;
      if (m_verbosity > 1) {
        std::cout << " Cutting track ID " << ((int)itrack)
          << "  because only has clusters in "
          << (in_0 ? "inner sector  (layers <23)"
            : in_1 ? "middle sector  (layers 23-39)"
            : "outer   (layers >40) ")
          << std::endl;
      }
      return true;
    }
  }
  return false;
}

void PHGhostRejection::find_ghosts(const std::vector<float>& trackChi2)
{
  if (m_verbosity > 0)
  {
    std::cout << "PHGhostRejection beginning track map size " << seeds.size() << std::endl;
  }

  // cut now pt tracks
  if (_min_pt>0.) {
    for (unsigned int i=0;i<seeds.size();++i) {
      if (seeds[i].get_pt() < _min_pt) {
        m_rejected[i] = true;
      }
    }
  }

  // Elimate low-interest track, and try to eliminate repeated tracks
  std::set<unsigned int> matches_set;
  std::multimap<unsigned int, unsigned int> matches;
  for (size_t trid1 = 0; trid1 < seeds.size(); ++trid1)
  {
    if (m_rejected[trid1]) { continue; }
    const auto& track1 = seeds[trid1];
    const float track1phi = track1.get_phi();

    const auto track1_pos = TrackSeedHelper::get_xyz(&track1);
    const float track1eta = track1.get_eta();
    for (size_t trid2 = trid1+1; trid2 < seeds.size(); ++trid2)
    {
      if (m_rejected[trid2])
      {
        continue;
      }

      const auto& track2 = seeds[trid2];
      const auto track2_pos = TrackSeedHelper::get_xyz(&track2);
      const float track2eta = track2.get_eta();
      auto delta_phi = std::abs(track1phi - track2.get_phi());

      if (delta_phi > 2 * M_PI) {
        delta_phi = delta_phi - 2*M_PI;
      }
      if (delta_phi < _phi_cut &&
          std::fabs(track1eta - track2eta) < _eta_cut &&
          std::fabs(track1_pos.x() - track2_pos.x()) < _x_cut &&
          std::fabs(track1_pos.y() - track2_pos.y()) < _y_cut &&
          std::fabs(track1_pos.z() - track2_pos.z()) < _z_cut)
      {
        matches_set.insert(trid1);
        matches.insert(std::pair(trid1, trid2));

        if (m_verbosity > 1)
        {
          std::cout << "Found match for tracks " << trid1 << " and " << trid2 << std::endl;
        }
      }
    }
  }

  for (auto set_it : matches_set)
  {
    if (m_rejected[set_it]) { continue; } // already rejected
    auto match_list = matches.equal_range(set_it);

    auto tr1 = seeds[set_it];
    double best_qual = trackChi2.at(set_it);
    unsigned int best_track = set_it;

    if (m_verbosity > 1)
    {
      std::cout << " ****** start checking track " << set_it << " with best quality " << best_qual << " best_track " << best_track << std::endl;
    }

    for (auto it = match_list.first; it != match_list.second; ++it)
    {
      if (m_verbosity > 1)
      {
        std::cout << "    match of track " << it->first << " to track " << it->second << std::endl;
      }

      auto tr2 = seeds[it->second];

      // Check that these two tracks actually share the same clusters, if not skip this pair
      bool is_same_track = checkClusterSharing(tr1, tr2);
      if (!is_same_track)
      {
        continue;
      }

      // which one has the best quality?
      double tr2_qual = trackChi2.at(it->second);
      if (m_verbosity > 1)
      {
        std::cout << "       Compare: best quality " << best_qual << " track 2 quality " << tr2_qual << std::endl;

        const auto track1_pos = TrackSeedHelper::get_xyz(&tr1);
        std::cout << "       tr1: phi " << tr1.get_phi() << " eta " << tr1.get_eta()
                  << " x " << track1_pos.x() << " y " << track1_pos.y() << " z " << track1_pos.z() << std::endl;

        const auto track2_pos = TrackSeedHelper::get_xyz(&tr2);
        std::cout << "       tr2: phi " << tr2.get_phi() << " eta " << tr2.get_eta()
                  << " x " << track2_pos.x() << " y " << track2_pos.y() << " z " << track2_pos.z() << std::endl;
      }

      if (tr2_qual < best_qual)
      {
        if (m_verbosity > 1)
        {
          std::cout << "       --------- Track " << it->second << " has better quality, erase track " << best_track << std::endl;
          std::cout << " rejecting track ID " << ((int)best_track) << "  because it is a ghost " << std::endl;
        }
        m_rejected[best_track] = true;
        best_qual = tr2_qual;
        best_track = it->second;
      }
      else
      {
        if (m_verbosity > 1)
        {
          std::cout << "       --------- Track " << best_track << " has better quality, erase track " << it->second << std::endl;
          std::cout << " rejecting track ID " << ((int)best_track) << "  because it is a ghost " << std::endl;
        }
        m_rejected[it->second] = true;
      }
    }
    if (m_verbosity > 1)
    {
      std::cout << " best track " << best_track << " best_qual " << best_qual << std::endl;
    }
  }

  // delete ghost tracks
  if (m_verbosity > 1)
  {
    for (unsigned int it=0 ; it<m_rejected.size(); ++it)
    {
      if (m_rejected[it]) {
        std::cout << " rejecting track ID " << it << std::endl;
      }
    }
  }

  if (m_verbosity > 0)
  {
    int n_ghost = std::count(m_rejected.begin(), m_rejected.end(), true);
    std::cout << " Track list sizes: n_init(" << ((int)m_rejected.size()) <<") - n_ghost("<<n_ghost << ") = n_good(" << (m_rejected.size()-n_ghost) << ")" << std::endl;
  }
}

// there is no check, at this point, about which is the best chi2 track
bool PHGhostRejection::checkClusterSharing(const TrackSeed& tr1, const TrackSeed& tr2) const
{
  // count shared clusters that tr1 and tr2 share many clusters
  size_t nclus_tr1 = tr1.size_cluster_keys();
  size_t nclus_tr2 = tr2.size_cluster_keys();
  size_t n_shared_clus = 0;

  for (auto key_tr1 = tr1.begin_cluster_keys();
       key_tr1 != tr1.end_cluster_keys();
       ++key_tr1)
  {
    if (tr2.find_cluster_key(*key_tr1) != tr2.end_cluster_keys())
    {
      ++n_shared_clus;
    }
  }

  if (m_verbosity > 2)
  {
    std::cout << " N-clusters tr1: " << nclus_tr1 << " N-clusters tr2: " << nclus_tr2 << " N-clusters shared: " << n_shared_clus << std::endl;
  }
  size_t nreq = 2 * n_shared_clus + 1;
  return (nreq > nclus_tr1) || (nreq > nclus_tr2);
}
