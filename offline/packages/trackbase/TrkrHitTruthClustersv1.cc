/**
 * @file trackbase/TrkrHitTruthClustersv1.cc
 * @author D. Stewart
 * @date April 2022
 * @brief Implementation of TrkrHitTruthClustersv1
 */

#include <math.h>
#include <iomanip>
#include "TrkrHitTruthClustersv1.h"

#include <algorithm>
#include <ostream>               // for operator<<, endl, basic_ostream, bas...

void
TrkrHitTruthClustersv1::Reset()
{ m_map.clear(); }

void
TrkrHitTruthClustersv1::print_clusters(std::ostream &os) const
{
  os << "-----TrkrHitTruthClustersv1-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;
  for( const auto& entry : m_map )
  {
    int track_id = entry.first;
    os << "   Centroids for track " << track_id << std::endl;
    os << std::setw(6);
    for (auto centroid : entry.second) {
      os << " layer: " << centroid.layer_id << " phi ("
         << std::setw(10) << centroid.phi_ave <<" pm "
         << std::setw(10) << centroid.phi_stdev<<") "
         << " z   ("
         << std::setw(10) << centroid.z_ave<<" pm "
         << std::setw(10) << centroid.z_stdev<<")  sum_E ("
         << std::setw(10) << centroid.sum_E << ") " << std::endl;
    }
    os << std::endl;
  }
  os << " end of embedded track energy centroid listings for this event" << std::endl;
  os << "---------------------------------------------------------------" << std::endl;
  return;
}

/* TrkrHitTruthClustersv1::CentroidsFor1Track& TrkrHitTruthClustersv1::add_track_centroids(const int i_track) { */
/*   //assumed that the tracks are added in increasing order */
/*   m_map[i_track] = {}; */
/*   return m_map[i_track]; */
/* } */
TrkrHitTruthClustersv1::VecEC& TrkrHitTruthClustersv1::get_new_centroids_vec(short track_id) {
  m_map[track_id] = {};
  return m_map[track_id];
}


