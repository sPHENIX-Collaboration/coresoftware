#include "TpcPolyClusterTrkrClusterConverter.h"

#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include "TpcCrossingDecision.h"
#include "TpcCrossingDecisionContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrDefs.h>

#include <tpc/TpcClusterMover.h>

#include <g4detectors/PHG4TpcGeomContainer.h>

#include <Acts/Definitions/Units.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace
{
  unsigned int adc_to_uint16(const double value)
  {
    if (!std::isfinite(value) || value <= 0.0) { return 0;
}
    return static_cast<unsigned int>(std::min(value, static_cast<double>(std::numeric_limits<unsigned short>::max())));
  }

  char size_to_char(const unsigned int value)
  {
    return static_cast<char>(std::min(value, 127U));
  }

  float finite_nonnegative_or_zero(const double value)
  {
    return std::isfinite(value) && value > 0.0 ? static_cast<float>(value) : 0.0F;
  }
}

TpcPolyClusterTrkrClusterConverter::TpcPolyClusterTrkrClusterConverter(const std::string& name)
  : SubsysReco(name)
{
}

TpcPolyClusterTrkrClusterConverter::~TpcPolyClusterTrkrClusterConverter() = default;

int TpcPolyClusterTrkrClusterConverter::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) { return Fun4AllReturnCodes::ABORTRUN;
}
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) { return Fun4AllReturnCodes::ABORTRUN;
}
  if (!initializeClusterMover(topNode)) { return Fun4AllReturnCodes::ABORTRUN;
}
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyClusterTrkrClusterConverter::getNodes(PHCompositeNode* topNode)
{
  m_polyClusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_polyClusterNodeName);
  if (!m_polyClusters)
  {
    std::cerr << Name() << "::getNodes - missing " << m_polyClusterNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_polyTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_polyTrackNodeName);
  if (!m_polyTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_polyTrackNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_crossingDecisions = findNode::getClass<TpcCrossingDecisionContainer>(topNode, m_crossingDecisionNodeName);
  if (!m_crossingDecisions)
  {
    std::cerr << Name() << "::getNodes - missing " << m_crossingDecisionNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_geometry)
  {
    std::cerr << Name() << "::getNodes - missing ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tpcGeomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!m_tpcGeomContainer)
  {
    std::cerr << Name() << "::getNodes - missing TPCGEOMCONTAINER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyClusterTrkrClusterConverter::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_outputClusters = findNode::getClass<TrkrClusterContainer>(topNode, m_outputNodeName);
  if (!m_outputClusters)
  {
    m_outputClusters = new TrkrClusterContainerv4();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_outputClusters, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TpcPolyClusterTrkrClusterConverter::clearOutputTpcClusters()
{
  if (!m_outputClusters) { return;
}

  const auto hitsetkeys = m_outputClusters->getHitSetKeys(TrkrDefs::TrkrId::tpcId);
  for (const auto hitsetkey : hitsetkeys)
  {
    m_outputClusters->removeClusters(hitsetkey);
  }
}

void TpcPolyClusterTrkrClusterConverter::buildTrackMap()
{
  m_tracksBySourceId.clear();
  if (!m_polyTracks) { return;
}

  for (unsigned int itrack = 0; itrack < m_polyTracks->size(); ++itrack)
  {
    const Tpc_PolyTrack* track = m_polyTracks->get_track(itrack);
    if (!isAcceptedTrack(track)) { continue;
}
    m_tracksBySourceId[track->get_source_assembled_track_id()] = track;
  }
}

bool TpcPolyClusterTrkrClusterConverter::isAcceptedTrack(const Tpc_PolyTrack* track) const
{
  if (!track || !track->isValid() || track->get_fit_status() == 0) { return false;
}
  if (track->get_nclusters() <= 20) { return false;
}

  const double pt = std::hypot(track->get_px(), track->get_py());
  return std::isfinite(pt) && pt > 0.1;
}

bool TpcPolyClusterTrkrClusterConverter::initializeClusterMover(PHCompositeNode* topNode)
{
  if (!m_geometry || !m_tpcGeomContainer || !topNode) { return false;
}

  m_clusterMover = std::make_unique<TpcClusterMover>();
  m_clusterMover->set_verbosity(Verbosity());
  m_clusterMover->initialize_geometry(m_tpcGeomContainer, m_geometry, topNode);
  return true;
}

TrkrDefs::cluskey TpcPolyClusterTrkrClusterConverter::getClusterKey(const Tpc_PolyCluster* cluster) const
{
  if (!cluster || cluster->size_hits() == 0) { return TrkrDefs::CLUSKEYMAX;
}

  const TrkrDefs::hitsetkey hitsetkey = cluster->get_hit_index(0).first;
  TrkrDefs::cluskey cluskey = cluster->get_trkr_cluster_key();
  const bool stored_key_is_tpc = cluskey != TrkrDefs::CLUSKEYMAX &&
                                 static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(cluskey)) == TrkrDefs::TrkrId::tpcId &&
                                 TrkrDefs::getHitSetKeyFromClusKey(cluskey) == hitsetkey;
  if (!stored_key_is_tpc)
  {
    cluskey = TrkrDefs::genClusKey(hitsetkey, cluster->get_cluster_id());
  }

  return cluskey;
}

bool TpcPolyClusterTrkrClusterConverter::localFromMovedGlobal(const Tpc_PolyCluster* cluster,
                                                              const TrkrDefs::cluskey cluskey,
                                                              const std::array<double, 3>& moved_global,
                                                              short crossing,
                                                              float& local_x,
                                                              float& local_y,
                                                              unsigned short& subsurfkey,
                                                              unsigned long long& surface_id) const
{
  if (!cluster || !m_geometry || cluster->size_hits() == 0) { return false;
}

  const TrkrDefs::hitsetkey hitsetkey = cluster->get_hit_index(0).first;
  const Acts::Vector3 global(moved_global[0], moved_global[1], moved_global[2]);
  if (!std::isfinite(global.x()) || !std::isfinite(global.y()) || !std::isfinite(global.z())) { return false;
}

  TrkrDefs::subsurfkey new_subsurfkey = 0;
  Surface surface = m_geometry->get_tpc_surface_from_coords(hitsetkey, global, new_subsurfkey);
  if (!surface)
  {
    const auto seed_iter = m_seedSubSurfKeys.find(cluskey);
    if (seed_iter == m_seedSubSurfKeys.end()) { return false;
}
    new_subsurfkey = seed_iter->second;
    surface = m_geometry->maps().getTpcSurface(hitsetkey, new_subsurfkey);
  }
  if (!surface) { return false;
}
  surface_id = surface->geometryId().value();

  const auto& geo_context = m_geometry->geometry().getGeoContext();
  Acts::Vector3 local = surface->localToGlobalTransform(geo_context).inverse() * (global * Acts::UnitConstants::cm);
  local /= Acts::UnitConstants::cm;
  const unsigned int side = TpcDefs::getSide(hitsetkey);
  const double drift_velocity = m_geometry->get_drift_velocity();
  if (drift_velocity <= 0.0 || !std::isfinite(drift_velocity)) { return false;
}

  const double half_drift = 0.5 * m_geometry->get_max_driftlength();
  const double z_bunch_separation = m_crossingPeriodNs * drift_velocity;
  const double zloc_uncorrected = (side == 0) ?
      local.y() + static_cast<double>(crossing) * z_bunch_separation :
      local.y() - static_cast<double>(crossing) * z_bunch_separation;
  const double tuncorrected = (side == 0) ? (zloc_uncorrected + half_drift) / drift_velocity : (half_drift - zloc_uncorrected) / drift_velocity;
  const double stored_t_uncorrected = tuncorrected - m_geometry->get_tpc_tzero() - m_geometry->get_sampa_tzero_bias();
  if (!std::isfinite(stored_t_uncorrected)) { return false;
}

  local_x = static_cast<float>(local.x());
  local_y = static_cast<float>(stored_t_uncorrected);
  subsurfkey = new_subsurfkey;
  return std::isfinite(local_x) && std::isfinite(local_y);
}

bool TpcPolyClusterTrkrClusterConverter::seedOutputCluster(const Tpc_PolyCluster* cluster, TrkrDefs::cluskey& cluskey)
{
  if (!cluster || !cluster->isValid() || !m_outputClusters) { return false;
}
  if (cluster->size_hits() == 0) { return false;
}

  const TrkrDefs::hitsetkey hitsetkey = cluster->get_hit_index(0).first;
  if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(hitsetkey)) != TrkrDefs::TrkrId::tpcId) { return false;
}

  cluskey = getClusterKey(cluster);
  if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(cluskey)) != TrkrDefs::TrkrId::tpcId) { return false;
}

  const Acts::Vector3 centroid(cluster->get_centroid_x(), cluster->get_centroid_y(), cluster->get_centroid_z());
  if (!std::isfinite(centroid.x()) || !std::isfinite(centroid.y()) || !std::isfinite(centroid.z())) { return false;
}

  TrkrDefs::subsurfkey subsurfkey = 0;
  Surface surface = m_geometry->get_tpc_surface_from_coords(hitsetkey, centroid, subsurfkey);
  if (!surface)
  {
    subsurfkey = 0;
  }
  m_seedSubSurfKeys[cluskey] = subsurfkey;

  auto out = std::make_unique<TrkrClusterv5>();
  out->setSubSurfKey(subsurfkey);
  out->setLocalX(0.0F);
  out->setLocalY(0.0F);
  m_outputClusters->addClusterSpecifyKey(cluskey, out.release());
  return true;
}

void TpcPolyClusterTrkrClusterConverter::buildMovedClusterMap()
{
  m_movedGlobals.clear();
  m_seedSubSurfKeys.clear();
  if (!m_polyClusters || !m_clusterMover) { return;
}

  std::map<unsigned int, std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>>> clusters_by_track;
  for (unsigned int icluster = 0; icluster < m_polyClusters->size(); ++icluster)
  {
    const Tpc_PolyCluster* cluster = m_polyClusters->get_cluster(icluster);
    if (!cluster || !cluster->isValid()) { continue;
}

    const auto track_iter = m_tracksBySourceId.find(cluster->get_source_assembled_track_id());
    if (track_iter == m_tracksBySourceId.end()) { continue;
}

    TrkrDefs::cluskey cluskey = TrkrDefs::CLUSKEYMAX;
    if (!seedOutputCluster(cluster, cluskey)) { continue;
}

    const Acts::Vector3 centroid(cluster->get_centroid_x(),
                                 cluster->get_centroid_y(),
                                 cluster->get_centroid_z());
    m_movedGlobals[cluskey] = {{centroid.x(), centroid.y(), centroid.z()}};
    clusters_by_track[track_iter->first].emplace_back(cluskey, centroid);
  }

  for (const auto& track_clusters : clusters_by_track)
  {
    const auto moved_globals = m_clusterMover->processTrack(track_clusters.second);
    for (const auto& [cluskey, moved_global] : moved_globals)
    {
      m_movedGlobals[cluskey] = {{moved_global.x(), moved_global.y(), moved_global.z()}};
    }
  }

  clearOutputTpcClusters();
}

bool TpcPolyClusterTrkrClusterConverter::publishCluster(const Tpc_PolyCluster* cluster,
                                                        const Tpc_PolyTrack* track) const
{
  if (!cluster || !cluster->isValid() || !track || !m_outputClusters) { return false;
}
  if (cluster->size_hits() == 0) { return false;
}

  const TrkrDefs::hitsetkey hitsetkey = cluster->get_hit_index(0).first;
  if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(hitsetkey)) != TrkrDefs::TrkrId::tpcId) { return false;
}

  const TrkrDefs::cluskey cluskey = getClusterKey(cluster);
  if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(cluskey)) != TrkrDefs::TrkrId::tpcId) { return false;
}

  const auto moved_iter = m_movedGlobals.find(cluskey);
  if (moved_iter == m_movedGlobals.end()) { return false;
}

  if (!m_crossingDecisions) { return false;
}
  const TpcCrossingDecision* crossing_decision =
      m_crossingDecisions->get_decision(track->get_source_assembled_track_id());
  if (!crossing_decision) { return false;
}
  const short crossing = crossing_decision->get_selected_crossing();

  float local_x = std::numeric_limits<float>::quiet_NaN();
  float local_y = std::numeric_limits<float>::quiet_NaN();
  unsigned short subsurfkey = 0;
  unsigned long long surface_id = 0;
  if (!localFromMovedGlobal(cluster, cluskey, moved_iter->second, crossing, local_x, local_y, subsurfkey, surface_id))
  {
    return false;
  }

  auto out = std::make_unique<TrkrClusterv5>();
  const unsigned int adc = adc_to_uint16(cluster->get_adc());
  out->setAdc(adc);
  out->setMaxAdc(static_cast<uint16_t>(adc));
  out->setEdge(0);
  out->setPhiSize(size_to_char(cluster->get_phi_width()));
  out->setZSize(size_to_char(cluster->get_time_width()));
  out->setSubSurfKey(subsurfkey);
  out->setOverlap(0);
  out->setLocalX(local_x);
  out->setLocalY(local_y);
  //out->setPhiError(0.002);
  //out->setPhiError(0.002);
  out->setPhiError(finite_nonnegative_or_zero(std::hypot(cluster->get_rms_x(), cluster->get_rms_y())));
  out->setZError(finite_nonnegative_or_zero(std::fabs(cluster->get_rms_z())));

  if (Verbosity() > 1)
  {
    const uint8_t sector = TpcDefs::getSectorId(hitsetkey);
    const uint8_t side = TpcDefs::getSide(hitsetkey);
    const uint8_t layer = TrkrDefs::getLayer(hitsetkey);
    const Acts::Vector3 trkr_global = m_geometry ? m_geometry->getGlobalPosition(cluskey, out.get()) : Acts::Vector3(
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::quiet_NaN());

    std::cout << Name() << "::publishCluster"
              << " track_id=" << track->get_track_id()
              << " source_track=" << track->get_source_assembled_track_id()
              << " crossing=" << crossing
              << " poly_cluster=" << cluster->get_cluster_id()
              << " hitsetkey=0x" << std::hex << hitsetkey
              << " cluskey=0x" << cluskey
              << " surface_id=0x" << surface_id << std::dec
              << " sector=" << static_cast<unsigned int>(sector)
              << " side=" << static_cast<unsigned int>(side)
              << " layer=" << static_cast<unsigned int>(layer)
              << " poly_pos=(" << cluster->get_centroid_x()
              << ", " << cluster->get_centroid_y()
              << ", " << cluster->get_centroid_z() << ")"
              << " moved_global=(" << moved_iter->second[0]
              << ", " << moved_iter->second[1]
              << ", " << moved_iter->second[2] << ")"
              << " trkr_local=(rphi,stored_t)=(" << out->getLocalX()
              << ", " << out->getLocalY() << ")"
              << " subsurfkey=" << static_cast<unsigned int>(out->getSubSurfKey())
              << " adc=" << out->getAdc()
              << " max_adc=" << out->getMaxAdc()
              << " phisize=" << out->getPhiSize()
              << " zsize=" << out->getZSize()
              << " edge=" << static_cast<unsigned int>(out->getEdge())
              << " overlap=" << static_cast<unsigned int>(out->getOverlap())
              << " phi_error=" << out->getRPhiError()
              << " z_error=" << out->getZError()
              << " trkr_global=(" << trkr_global.x()
              << ", " << trkr_global.y()
              << ", " << trkr_global.z() << ")"
              << " dz_check=" << trkr_global.z() - moved_iter->second[2]
              << std::endl;
  }

  m_outputClusters->addClusterSpecifyKey(cluskey, out.release());
  return true;
}

int TpcPolyClusterTrkrClusterConverter::process_event(PHCompositeNode* topNode)
{
  if (!m_polyClusters || !m_polyTracks || !m_outputClusters || !m_crossingDecisions || !m_geometry || !m_tpcGeomContainer || !m_clusterMover)
  {
    if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        !initializeClusterMover(topNode))
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  clearOutputTpcClusters();
  buildTrackMap();
  buildMovedClusterMap();

  unsigned int nconverted = 0;
  unsigned int nmissingTrack = 0;
  unsigned int nfailed = 0;
  for (unsigned int icluster = 0; icluster < m_polyClusters->size(); ++icluster)
  {
    const Tpc_PolyCluster* cluster = m_polyClusters->get_cluster(icluster);
    if (!cluster || !cluster->isValid()) { continue;
}

    const auto track_iter = m_tracksBySourceId.find(cluster->get_source_assembled_track_id());
    if (track_iter == m_tracksBySourceId.end())
    {
      ++nmissingTrack;
      continue;
    }

    if (publishCluster(cluster, track_iter->second)) { ++nconverted;
    } else { ++nfailed;
}
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - poly_clusters=" << m_polyClusters->size()
              << " converted=" << nconverted
              << " missing_track=" << nmissingTrack
              << " failed=" << nfailed << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
