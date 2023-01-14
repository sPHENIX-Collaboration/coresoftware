/**
 * @file trackbase/TrkrTruthClustersv1.cc
 * @author D. Stewart
 * @date April 2022
 * @brief Implementation of TrkrTruthClustersv1
 */
#include "TrkrTruthClustersv1.h"
#include "TrkrClusterv4.h"

#include <trackbase/TrkrDefs.h>

#include <TString.h>

#include <array>
#include <iomanip>
#include <set>
#include <algorithm>


TrkrTruthClusters::ConstRange TrkrTruthClustersv1::getClusters (short trkId) const {
  if (trkId==-1) {
    return { m_data.begin(), m_data.end() };
  }
  auto begin = std::lower_bound(m_data.begin(), m_data.end(), trkId, CompTrkId());
  auto end   = std::upper_bound(m_data.begin(), m_data.end(), trkId, CompTrkId());
  return { begin, end };
}

std::vector<short> TrkrTruthClustersv1::getTrkIds (short layer) const {
  // get all track ids with a given layer
  std::vector<short> rvec {};
  auto iter = m_data.begin();
  auto all_end = m_data.end();
  while (iter != all_end) {
    short trkid = iter->first;
    auto this_end = std::upper_bound(iter,all_end,trkid,CompTrkId());
    if (std::binary_search(iter,this_end,layer,CompLayer())) { 
      rvec.push_back(trkid); 
    }
    iter = this_end;
  }
  return rvec;
}

std::vector<short> TrkrTruthClustersv1::getLayers (short trkid) const {
  std::set<short> layers{};
  auto range = getClusters(trkid); // will loop over one track or (by default) all tracks
  for (auto& E = range.first; E != range.second; ++E) {
    layers.insert(TrkrDefs::getLayer(E->second));
  }
  std::vector<short> rvec;
  for (auto& layer : layers) rvec.push_back(layer);
  return rvec;
}

bool TrkrTruthClustersv1::hasTrkId (short trkid) const {
  return std::binary_search(m_data.begin(), m_data.end(), trkid, CompTrkId());
}

// must simply search through all entries...
bool TrkrTruthClustersv1::hasTrkIdLayer (short trkId, short layer) const {
  return std::binary_search(m_data.begin(), m_data.end(), 
      std::pair<short,short>{trkId,layer}, CompTrkIdLayer());
}

// must search through all entries...
// if this must be done many times, then user getLayerss once, and then search that vector
bool TrkrTruthClustersv1::hasLayer (short layer) const {
  auto iter = m_data.begin();
  auto end_all = m_data.end();
  while (iter != end_all) {
    short trkid = iter->first;
    auto end_loc = std::upper_bound(iter,end_all,trkid,CompTrkId());
    if (std::binary_search(iter, end_loc, layer, CompLayer()))
    { return true; }
    iter = end_loc;
  }
  return false;
}

void TrkrTruthClustersv1::addTruthCluster (short trkid, TrkrDefs::cluskey clusterkey) {
  m_data.push_back({trkid,clusterkey});
}

//void TrkrTruthClustersv1::identify(std::ostream &os) const {
  //os << " not null " << std::endl;
//}

//void TrkrTruthClustersv1::identify(std::ostream &os, TrkrClusterContainer* container) const {
  //if (container) os << " not null " << std::endl;
//}

std::ostream& 
TrkrTruthClustersv1::print_clusters(TrkrClusterContainer* clusters, std::ostream &os) const
{
  os << "-----TrkrTruthClustersv1-----" << std::endl;

  os << "Tracks ids present: " << m_data.size() << std::endl;
  for (auto trk : getTrkIds()) os << "  " << trk << std::endl;

  os << "Layer ids present: " << m_data.size() << std::endl;
  auto layers = getLayers();
  for (auto layer : getLayers()) os << "  " << layer << std::endl;
  os << Form(" %10s %8s %14s %10s %12s %12s %12s %12s","TrkrId","TpcLayer",
      "hit data:","neff_electrons","phi","Z","phibin_range","Zbin_range") << std::endl;

  os << " Track Layer's present (v1): " << std::endl;
  auto range = getClusters();
  for (auto entry = range.first; entry != range.second; ++entry) {
    os << "Track id: " << entry->first << " layer: " << TrkrDefs::getLayer(entry->second) << std::endl;
    auto cluster = clusters->findCluster(entry->second);
    os << Form(" %10i %8i %14s %10i %12f %12f %12f %12f", entry->first, TrkrDefs::getLayer(entry->second), "", 
        cluster->getAdc(), cluster->getPosition(0), cluster->getPosition(1), 
        cluster->getPhiSize(), cluster->getZSize()) << std::endl;
    os << " is present? track(" << hasTrkId(entry->first) 
      << " hasLayers(" << hasLayer(TrkrDefs::getLayer(entry->second))
      << " hasBoth("    << hasTrkIdLayer(entry->first,TrkrDefs::getLayer(entry->second)) << std::endl;
  }

  os << " Track Layer's present (v2): " << std::endl;
  auto trackids = getTrkIds();
  os << Form(" %10s %8s %14s %10s %12s %12s %12s %12s","TrkrId","TpcLayer",
      "hit data:","neff_electrons","phi","Z","phibin_range","Zbin_range") << std::endl;
  for (auto trackid : trackids) { //range.first; entry != range.second; ++entry) {
    os << "Track id: " << trackid;
    os << " -> Layer ids: " << std::endl; 
    for (auto layer : getLayers(trackid)) {
      os << "  " << layer << std::endl;

      auto clusid = getCluskey(trackid, layer);
      auto cluster = clusters->findCluster(clusid);
      os << Form(" %10i %8i %14s %10i %12f %12f %12f %12f", trackid, layer, "", 
          cluster->getAdc(), cluster->getPosition(0), cluster->getPosition(1), 
          cluster->getPhiSize(), cluster->getZSize()) << std::endl;
    }
    /* os << " is present? track(" << hasTrkId(entry->first) << */ 
      /* << " hasLayers(" << hasLayer(TrkrDefs::getLayer(entry->second) */ 
      /* << " hasBoth("    << hasTrkIdLayer(entry->first,TrkrDefs::getLayer(entry->second)) << std::endl; */
  }

  return os;
}

TrkrDefs::cluskey TrkrTruthClustersv1::getCluskey(short trkId, short layer) const {
  auto iter = lower_bound(m_data.begin(), m_data.end(), 
      std::pair<short,short>{trkId,layer}, CompTrkIdLayer());
  if (iter == m_data.end()) 
  { return 0; }

  return iter->second;
}
