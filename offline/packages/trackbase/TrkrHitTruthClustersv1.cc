/**
 * @file trackbase/TrkrHitTruthClustersv1.cc
 * @author D. Stewart
 * @date April 2022
 * @brief Implementation of TrkrHitTruthClustersv1
 */
#include "TrkrHitTruthClustersv1.h"
#include "TrkrClusterv4.h"

#include <iomanip>
#include <set>
#include <algorithm>
#include "TString.h"


TrkrHitTruthClusters::ConstRange TrkrHitTruthClustersv1::getClusters (short trkId) const {
  if (trkId==-1) {
  /* if (false) std::cout << trkId << std::endl; */ 
    return { m_data.begin(), m_data.end() };
  }
  auto begin = std::lower_bound(m_data.begin(), m_data.end(), trkId, 
      [](const Entry& lhs, const short& rhs){ return lhs.first.first < rhs; });
  if (begin == m_data.end() || begin->first.first != trkId) return { m_data.end(), m_data.end() };
  auto end   = std::lower_bound(m_data.begin(), m_data.end(), trkId+1, 
      [](const Entry& lhs, const short& rhs){ return lhs.first.first < rhs; });
  return { begin, end };
}

std::vector<short> TrkrHitTruthClustersv1::getTrkIds (short layer) const {
  std::vector<short> rvec {};
  if (layer == -1) {
    short trkid =-1;
    for (auto& E : m_data) {
      if (E.first.first == trkid) {
        trkid = E.first.first;
        rvec.push_back(trkid);
      }
    }
  } else {
    for (auto& E : m_data) {
      if (E.first.second == layer) {
        rvec.push_back(E.first.first);
      }
    }
  }
  return rvec;
}

std::vector<short> TrkrHitTruthClustersv1::getLayerIds (short trkid) const {
  std::set<short> layers{};
  auto range = getClusters(trkid); // will loop over one track or (by default) all tracks
  for (auto& E = range.first; E != range.second; ++E) {
    layers.insert(E->first.second);
  }

  std::vector<short> rvec;
  for (auto& layer : layers) rvec.push_back(layer);
  return rvec;
}

bool TrkrHitTruthClustersv1::hasTrkId (short trkId) const {
  auto begin = std::lower_bound(m_data.begin(), m_data.end(), trkId, 
      [](const Entry& lhs, const short& rhs){ return lhs.first.first < rhs; });
  return (begin != m_data.end() && begin->first.first == trkId);
}

// must simply search through all entries...
bool TrkrHitTruthClustersv1::hasTrkIdLayerId (short trkId, short layerId) const {
  Key key {trkId, layerId}; 
  return std::lower_bound(m_data.begin(), m_data.end(), key, [](const Entry& lhs, const Key& rhs) 
      { return lhs.first<rhs;}) == m_data.end();
}

// must search through all entries...
// if this must be done many times, then user getLayerIds once, and then search that vector
bool TrkrHitTruthClustersv1::hasLayerId (short layerid) const {
  for (auto& E : m_data) {
    if (E.first.second == layerid) return true;
  }
  return false;
}

void TrkrHitTruthClustersv1::addTruthCluster (short trkid, MapToPadPlanePassData& _hit) {
  if (_hit.layer == 0) return; // there is no-zero'th layer; do nothing
  Key key { trkid, _hit.layer };
  TrkrClusterv4* cluster = new TrkrClusterv4();
  cluster->setPosition( 0, _hit.phi_integral  / _hit.neff_electrons);
  cluster->setPosition( 1, _hit.time_integral / _hit.neff_electrons);
  cluster->setPhiSize ( _hit.phi_bin_hi  - _hit.phi_bin_lo  +1);
  cluster->setZSize   ( _hit.time_bin_hi - _hit.time_bin_lo +1);
  cluster->setAdc     ( _hit.neff_electrons);
  m_data.push_back ( {key, cluster} );
  _hit.reset();
}

void TrkrHitTruthClustersv1::Reset() {
  for (auto& E : m_data) delete E.second;
  m_data.clear();
}

void
TrkrHitTruthClustersv1::print_clusters(std::ostream &os) const
{
  os << "-----TrkrHitTruthClustersv1-----" << std::endl;
  os << "Number of associations: " << m_data.size() << std::endl;
  // TrkId TpcLayer  Cluster: neff_electrons phi Z phibin_range Zbin_range" 
  os << Form(" %10s %8s %14s %10s %12s %12s %12s %12s","TrkrId","TpcLayer",
      "hit data:","neff_electrons","phi","Z","phibin_range","Zbin_range") << std::endl;
  for (auto& E : m_data) {
    auto key = E.first;
    auto hit = E.second;
    os << Form(" %10i %8i %14s %10i %12f %12f %12f %12f", key.first, key.second, "", 
        hit->getAdc(), hit->getPosition(0), hit->getPosition(1), hit->getPhiSize(), hit->getZSize()) << std::endl;
  }
}

