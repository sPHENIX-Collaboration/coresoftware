/**
 * @file trackbase/TrkrHitTruthClustersv1.cc
 * @author D. Stewart
 * @date April 2022
 * @brief Implementation of TrkrHitTruthClustersv1
 */

#include "TrkrHitTruthClustersv1.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>  // for operator<<, endl, basic_ostream, bas...

void TrkrHitTruthClustersv1::Reset()
{
  m_map.clear();
}

void TrkrHitTruthClustersv1::print_clusters(std::ostream& os) const
{
  os << "-----TrkrHitTruthClustersv1-----" << std::endl;
  os << "Number of associations: " << m_map.size() << std::endl;

  for (const auto& entry : m_map)
  {
    int track_id = entry.first;
    os << "   Hits for track " << track_id << std::endl;
    os << std::setw(3) << "#" << " ";
    for (auto w : std::vector<const char*>{
             "phi_mean", "phi_stdev", "R_mean", "R_stddev",
             "z_mean", "z_std"})
    {
      os << " " << std::setw(10) << w;
    }
    os << std::endl;
    for (unsigned int i_cluster = 0; i_cluster < entry.second.size(); ++i_cluster)
    {
      os << std::setw(3) << i_cluster << " ";
      for (int k{0}; k < 6; ++k)
      {
        os << " " << std::setw(10) << entry.second[i_cluster][k];
      }
      os << std::endl;
    }
    os << std::endl;
  }
  os << "------------------------------" << std::endl;
  return;
}

void TrkrHitTruthClustersv1::push_truth_cluster(const int i_track,
                                                const std::array<double, 8>& phiRz_data, const double sum_E)
{
  /* phiRz_data (the input) is:
     *   sum(phi), sum(phi^2) 
     *   sum(R)  , sum(R^2)
     *   sum(eta), sum(eta^2)
     *   sum(phi_rotate), sum(phi_rotate^2)
     * The purpose of phi_rotate is to account for the cirular boundary conditions,
     * it is summed with each entry rotated by +pi if the input is negative
     */

  bool print_local = false;
  if (sum_E == 0)
  {
    if (print_local)
    {
      std::cout << "Trying to push TrkrHitTruthClusterv1 for embedded track \""
                << i_track << "\" with 0 entries" << std::endl
                << " -> Skipping entry" << std::endl;
    }
    return;
  }

  if (print_local)
  {
    std::cout << " pushing truth cluster: " << std::endl;
    for (auto& v : phiRz_data)
    {
      std::cout << " " << v;
    }
    std::cout << " sum_E: " << sum_E << std::endl;
  }
  std::array<double, 8> phiRz_calc{};
  for (int i : std::vector<int>{0, 2, 4, 6})
  {
    double Hj = phiRz_data[i];
    double Ej = phiRz_data[i + 1];
    double Wj = sum_E;
    double hj = Hj / Wj;
    phiRz_calc[i] = hj;                                // energy weighted mean
    phiRz_calc[i + 1] = std::sqrt(Ej / Wj - hj * hj);  // energy weighted std.
    if (print_local)
    {
      std::cout << " Added: (" << i << ") : " << phiRz_calc[i] << " " << phiRz_calc[i + 1] << std::endl;
    }
  }
  std::array<float, 6> phiRz_entry{};
  for (int i : std::vector<int>{2, 3, 4, 5})
  {
    phiRz_entry[i] = phiRz_calc[i];
  }

  double phi_mean{phiRz_calc[0]};
  double phi_std{phiRz_calc[1]};
  if (phi_std > phiRz_calc[7])
  {
    phi_std = phiRz_calc[7];
    phi_mean = phiRz_calc[6];
    if (phi_mean > M_PI)
    {
      phi_mean -= 2 * M_PI;
    }
  }
  phiRz_entry[0] = phi_mean;
  phiRz_entry[1] = phi_std;

  if (m_map.find(i_track) != m_map.end())
  {
    m_map[i_track].push_back(phiRz_entry);
  }
  else
  {
    m_map[i_track] = {phiRz_entry};
  }

  if (print_local)
  {
    for (auto& v_pairs : m_map)
    {
      std::cout << " track " << v_pairs.first << "  n-clusters: " << v_pairs.second.size() << std::endl;
    }
  }
}
