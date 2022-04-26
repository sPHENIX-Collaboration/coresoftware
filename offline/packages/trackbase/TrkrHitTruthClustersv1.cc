/**
 * @file trackbase/TrkrHitTruthClustersv1.cc
 * @author D. Stewart
 * @date April 2022
 * @brief Implementation of TrkrHitTruthClustersv1
 */

#include <math.h>
#include <iomanip>
#include "TrkrHitTruthClustersv1.h"
/* #include "TrkrDefs.h" */

/* #include <g4main/PHG4HitDefs.h> */

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

  // save the state of os before formating:
  //  (It would be nice to use ROOT's `Format` function, but lacking that...)
  std::ios old_state(nullptr);
  old_state.copyfmt(os); 
  os << std::setw(9); // for the purposes of the output statements


  for( const auto& entry : m_map )
  {
    int track_id = entry.first;
    os << "   Hits for track " << track_id << std::endl;
    os << "phi_mean" << " " << "phi_stddev" << " "
       << "eta_mean" << " " << "eta_stddev" << " " 
       << "z_mean"   << " " << "z_stddev" << " " << std::endl;
    for (unsigned int i_cluster=0; i_cluster<entry.second.size(); ++i_cluster) {
        os << i_cluster;
        for (int k{0}; k<6; ++k) os << entry.second[i_cluster][k];
        os << std::endl;
    }
    os << std::endl << std::endl;
  }
  os << "------------------------------" << std::endl;
  os.copyfmt(old_state);

  return;
}

void TrkrHitTruthClustersv1::push_truth_cluster(const int i_track, 
        const std::array<double,6>& phi_eta_z_data, const double sum_E)
{
    // polar map has phi, R, Z coordinates (mean and std for each) in an array
    std::array<float,6> phi_eta_z_entry;
    for (int i : std::vector<int>{0,2,4}) {
        double Hj = phi_eta_z_data[i];
        double Ej = phi_eta_z_data[i+1];
        double Wj = sum_E;
        phi_eta_z_entry[i]   = Hj/Wj; // energy weighted mean
        phi_eta_z_entry[i+1] = std::sqrt(Ej/Wj-Hj*Hj); // energy weighted std.
    }
    if (m_map.find(i_track) != m_map.end()) m_map[i_track].push_back(phi_eta_z_entry);
    else m_map[i_track] = {phi_eta_z_entry};
}
