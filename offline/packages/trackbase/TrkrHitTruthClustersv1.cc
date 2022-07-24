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
    os << std::setw(3) << "#" << " ";
    for (int i=0;i<55;++i) {
      const EnergyCentroid& centroid = entry.second[i];
      /* if (centroid.sum_E == 0) continue; */
      os << " layer: " << (i+1) << " phi ("<<centroid.phi_ave<<" pm "<<centroid.phi_stdev<<") "
                                       << " z   ("<<centroid.z_ave<<" pm "<<centroid.z_stdev<<")  sum_E ("<<
                                       centroid.sum_E << ") " << std::endl;
    }
    os << std::endl;
  }
  os << "------------------------------" << std::endl;
  return;
}

TrkrHitTruthClustersv1::CentroidsFor1Track& TrkrHitTruthClustersv1::add_track_centroids(const int i_track) {
  //assumed that the tracks are added in increasing order
  m_map[i_track] = {};
  return m_map[i_track];
}
