#include "TpcCrossingFinder.h"

#include "IdealPadMap.h"
#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackContainer.h"
#include "TpcCrossingDecisionContainerv1.h"
#include "TpcCrossingDecisionv1.h"
#include "Tpc_FittingTools.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <ffamodules/CDBInterface.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>
#include <trackbase/InttDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <phgarfield/PHGarfield.h>
#include <TPolyLine3D.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace
{
  constexpr unsigned int FirstLayer = 7;
  constexpr unsigned int LastLayer = 54;
  constexpr unsigned int NLayers = LastLayer - FirstLayer + 1;
  constexpr unsigned int NSides = 2;
  constexpr unsigned int NSectors = 12;
  constexpr double PhiConsistencyTolerance = 1.0e-10;
  constexpr unsigned char InvalidConfidenceTier = std::numeric_limits<unsigned char>::max();

  double wrap_phi(double phi)
  {
    while (phi > M_PI) { phi -= 2.0 * M_PI;
}
    while (phi <= -M_PI) { phi += 2.0 * M_PI;
}
    return phi;
  }

  double unwrap_phi_near(double phi, const double reference)
  {
    while (phi - reference > M_PI) { phi -= 2.0 * M_PI;
}
    while (phi - reference <= -M_PI) { phi += 2.0 * M_PI;
}
    return phi;
  }

  double clamp_unit(const double value)
  {
    return std::max(0.0, std::min(1.0, value));
  }

  double phi_sample_fraction(const unsigned int sample)
  {
    if (TpcCrossingFinder::NPhiSamples <= 1U) { return 0.0;
}
    return static_cast<double>(sample) / static_cast<double>(TpcCrossingFinder::NPhiSamples - 1U);
  }

  float finite_float_or_nan(const double value)
  {
    return std::isfinite(value) ? static_cast<float>(value) : std::numeric_limits<float>::quiet_NaN();
  }
}

TpcCrossingFinder::TpcCrossingFinder(const std::string& name)
  : SubsysReco(name)
{
}

TpcCrossingFinder::~TpcCrossingFinder()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
  delete m_garfield;
  m_garfield = nullptr;
}

int TpcCrossingFinder::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) { return Fun4AllReturnCodes::ABORTRUN;
}
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) { return Fun4AllReturnCodes::ABORTRUN;
}

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - failed to load IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4TpcGeom* layergeom = m_geomContainerTpc->GetLayerCellGeom(20);
  if (layergeom)
  {
    const double rot_x = layergeom->get_rot_x();
    const double rot_y = layergeom->get_rot_y();
    const double rot_z = layergeom->get_rot_z();
    const double place_x = layergeom->get_place_x();
    const double place_y = layergeom->get_place_y();
    const double place_z = layergeom->get_place_z();
    if (use_survey_geometry)
    {
      m_tpcMove = {{place_x, place_y, place_z}};
      m_tpcRotations = {{{{rot_x, rot_y, rot_z}}, {{0.0, 0.0, 0.0}}}};
    }
  }

  delete m_garfield;
  const std::string electricFieldMap = CDBInterface::instance()->getUrl("Tpc_PolySeeding_Efield");
  m_garfield = new PHGarfield(Name() + "_PHGarfield", electricFieldMap, m_kEffSide0, m_kEffSide1);
  configure_garfield(m_garfield);
  if (m_garfield->InitRun(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << Name() << "::InitRun - PHGarfield InitRun failed" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!build_drift_lookup()) { return Fun4AllReturnCodes::ABORTRUN;
}

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCrossingFinder::getNodes(PHCompositeNode* topNode)
{
  m_assembledTracks = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_inputNodeName);
  if (!m_assembledTracks)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterMap && Verbosity() > 0)
  {
    std::cout << Name() << "::getNodes - optional TRKR_CLUSTER node not found" << std::endl;
  }

  m_geomContainerTpc = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!m_geomContainerTpc)
  {
    std::cerr << Name() << "::getNodes - missing TPCGEOMCONTAINER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexMapNodeName);
  if (!m_vertexMap && Verbosity() > 0)
  {
    std::cout << Name() << "::getNodes - optional " << m_vertexMapNodeName << " node not found" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcCrossingFinder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_decisions = findNode::getClass<TpcCrossingDecisionContainerv1>(topNode, m_outputNodeName);
  if (!m_decisions)
  {
    m_decisions = new TpcCrossingDecisionContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_decisions, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void TpcCrossingFinder::configure_garfield(PHGarfield* garfield) const
{
  if (!garfield) { return;
}

  garfield->MoveTpc(m_tpcMove[0], m_tpcMove[1], m_tpcMove[2]);
  for (const auto& rotation : m_tpcRotations)
  {
    garfield->RotateTpc(rotation[0], rotation[1], rotation[2]);
  }
  garfield->SetCMVoltageDefault(m_cmVoltageDefault);
}

unsigned int TpcCrossingFinder::drift_lookup_index(const unsigned int layer_index,
                                                   const unsigned int side,
                                                   const unsigned int sector,
                                                   const unsigned int sample)
{
  return (((layer_index * NSides + side) * NSectors + sector) * NPhiSamples + sample);
}

bool TpcCrossingFinder::build_drift_lookup()
{
  if (!m_idealPadMap || !m_garfield) { return false;
}
  if (m_reverseDriftStepNs <= 0.0 || !std::isfinite(m_reverseDriftStepNs)) { return false;
}

  m_maxLookupTimeNs = std::numeric_limits<double>::max();
  for (DriftPolyline& polyline : m_driftLookup)
  {
    polyline.phi = 0.0;
    polyline.points.clear();
  }

  unsigned int nbuilt = 0;
  for (unsigned int layer = FirstLayer; layer <= LastLayer; ++layer)
  {
    const unsigned int layer_index = layer - FirstLayer;
    const double radius = m_idealPadMap->get_radius(layer);
    const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(layer);
    if (!std::isfinite(radius) || pads_per_sector == 0U) { return false;
}

    for (unsigned int side = 0; side < NSides; ++side)
    {
      const double z0 = (side == 0U) ? m_startZSouth : m_startZNorth;
      for (unsigned int sector = 0; sector < NSectors; ++sector)
      {
        for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
        {
          const unsigned int local_phibin = static_cast<unsigned int>(std::llround(
              phi_sample_fraction(sample) * static_cast<double>(pads_per_sector - 1U)));
          const unsigned int global_pad = sector * pads_per_sector + local_phibin;
          const double phi_local = m_idealPadMap->get_phi(side, sector, layer, local_phibin);
          const double phi_global = m_idealPadMap->get_phi(side, layer, global_pad);
          if (!std::isfinite(phi_local) || !std::isfinite(phi_global)) { return false;
}
          if (std::abs(wrap_phi(phi_global - phi_local)) > PhiConsistencyTolerance) { return false;
}

          const double x0 = radius * std::cos(phi_local);
          const double y0 = radius * std::sin(phi_local);
          TPolyLine3D* drift = m_garfield->ReverseDrift(x0, y0, z0, m_reverseDriftStepNs);
          if (!drift || drift->GetN() <= 0)
          {
            delete drift;
            return false;
          }

          const int npoints = drift->GetN();
          const Float_t* xyz = drift->GetP();
          if (!xyz || npoints <= 0)
          {
            delete drift;
            return false;
          }

          DriftPolyline& polyline = m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
          polyline.phi = phi_local;
          polyline.points.resize(static_cast<std::size_t>(npoints));
          for (int ipoint = 0; ipoint < npoints; ++ipoint)
          {
            const int idx = 3 * ipoint;
            DriftPoint& point = polyline.points[static_cast<std::size_t>(ipoint)];
            const double output_r = std::hypot(static_cast<double>(xyz[idx]), static_cast<double>(xyz[idx + 1]));
            const double output_phi = unwrap_phi_near(std::atan2(static_cast<double>(xyz[idx + 1]), static_cast<double>(xyz[idx])), phi_local);
            point.delta_r = static_cast<float>(output_r - radius);
            point.delta_phi = static_cast<float>(output_phi - phi_local);
            point.z = xyz[idx + 2];
          }
          m_maxLookupTimeNs = std::min(m_maxLookupTimeNs, static_cast<double>(npoints - 1) * m_reverseDriftStepNs);
          delete drift;
          ++nbuilt;
        }
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::build_drift_lookup - built " << nbuilt
              << " drift polylines max_lookup_time_ns=" << m_maxLookupTimeNs << std::endl;
  }
  return nbuilt == NLayers * NSides * NSectors * NPhiSamples && std::isfinite(m_maxLookupTimeNs) && m_maxLookupTimeNs > 0.0;
}

bool TpcCrossingFinder::sample_drift_lookup(const unsigned int layer,
                                            const unsigned int side,
                                            const unsigned int pad,
                                            const unsigned int tbin,
                                            const short crossing,
                                            double& x,
                                            double& y,
                                            double& z) const
{
  if (!m_idealPadMap) { return false;
}
  if (layer < FirstLayer || layer > LastLayer) { return false;
}
  if (side >= NSides) { return false;
}
  if (m_reverseDriftStepNs <= 0.0 || !std::isfinite(m_reverseDriftStepNs)) { return false;
}

  const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(layer);
  if (pads_per_sector == 0U) { return false;
}
  const unsigned int sector = pad / pads_per_sector;
  if (sector >= NSectors) { return false;
}

  const double hit_radius = m_idealPadMap->get_radius(layer);
  const double hit_phi = m_idealPadMap->get_phi(side, layer, pad);
  if (!std::isfinite(hit_radius) || !std::isfinite(hit_phi)) { return false;
}

  const double target_time_ns = (static_cast<double>(tbin) - m_t0) * m_tpcAdcClock
    - static_cast<double>(crossing) * m_crossingPeriodNs;
  if (target_time_ns <= 0.0 || !std::isfinite(target_time_ns)) { return false;
}

  const unsigned int layer_index = layer - FirstLayer;
  std::array<const DriftPolyline*, NPhiSamples> samples{};
  std::array<double, NPhiSamples> sample_phi{};
  for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
  {
    samples[sample] = &m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
    if (samples[sample]->points.empty()) { return false;
}
    sample_phi[sample] = samples[sample]->phi;
    if (sample > 0U) { sample_phi[sample] = unwrap_phi_near(sample_phi[sample], sample_phi[sample - 1U]);
}
  }

  const double unwrapped_hit_phi = unwrap_phi_near(hit_phi, sample_phi[NPhiSamples / 2U]);
  const bool increasing = sample_phi[NPhiSamples - 1U] >= sample_phi[0];

  bool bracket_found = false;
  unsigned int sample0 = 0;
  unsigned int sample1 = 0;
  double phi_fraction = 0.0;
  if ((increasing && unwrapped_hit_phi <= sample_phi[0]) || (!increasing && unwrapped_hit_phi >= sample_phi[0]))
  {
    bracket_found = true;
  }
  else if ((increasing && unwrapped_hit_phi >= sample_phi[NPhiSamples - 1U]) ||
           (!increasing && unwrapped_hit_phi <= sample_phi[NPhiSamples - 1U]))
  {
    sample0 = NPhiSamples - 1U;
    sample1 = NPhiSamples - 1U;
    bracket_found = true;
  }
  else
  {
    for (unsigned int sample = 0; sample + 1U < NPhiSamples; ++sample)
    {
      const bool in_interval = increasing ?
        (unwrapped_hit_phi >= sample_phi[sample] && unwrapped_hit_phi <= sample_phi[sample + 1U]) :
        (unwrapped_hit_phi <= sample_phi[sample] && unwrapped_hit_phi >= sample_phi[sample + 1U]);
      if (!in_interval) { continue;
}

      sample0 = sample;
      sample1 = sample + 1U;
      const double denom = sample_phi[sample1] - sample_phi[sample0];
      phi_fraction = (denom != 0.0) ? clamp_unit((unwrapped_hit_phi - sample_phi[sample0]) / denom) : 0.0;
      bracket_found = true;
      break;
    }
  }
  if (!bracket_found) { return false;
}

  auto sample_time = [this, target_time_ns](const DriftPolyline& polyline,
                                            double& delta_r,
                                            double& delta_phi,
                                            double& point_z) -> bool
  {
    const int npoints = static_cast<int>(polyline.points.size());
    if (npoints <= 0) { return false;
}
    const double max_time_ns = static_cast<double>(npoints - 1) * m_reverseDriftStepNs;
    if (target_time_ns > max_time_ns) { return false;
}

    const double fbin = target_time_ns / m_reverseDriftStepNs;
    const int i0 = std::min(static_cast<int>(std::floor(fbin)), npoints - 1);
    const int i1 = std::min(i0 + 1, npoints - 1);
    const double frac = fbin - static_cast<double>(i0);

    const DriftPoint& p0 = polyline.points[static_cast<std::size_t>(i0)];
    const DriftPoint& p1 = polyline.points[static_cast<std::size_t>(i1)];
    const double dphi0 = static_cast<double>(p0.delta_phi);
    const double dphi1 = unwrap_phi_near(static_cast<double>(p1.delta_phi), dphi0);
    delta_r = static_cast<double>(p0.delta_r) + frac * static_cast<double>(p1.delta_r - p0.delta_r);
    delta_phi = dphi0 + frac * (dphi1 - dphi0);
    point_z = static_cast<double>(p0.z) + frac * static_cast<double>(p1.z - p0.z);
    return std::isfinite(delta_r) && std::isfinite(delta_phi) && std::isfinite(point_z);
  };

  double delta_r0 = 0.0;
  double delta_phi0 = 0.0;
  double z0 = 0.0;
  const bool valid0 = sample_time(*samples[sample0], delta_r0, delta_phi0, z0);
  double delta_r1 = 0.0;
  double delta_phi1 = 0.0;
  double z1 = 0.0;
  const bool same_sample = sample0 == sample1;
  const bool valid1 = same_sample ? valid0 : sample_time(*samples[sample1], delta_r1, delta_phi1, z1);
  if (same_sample)
  {
    delta_r1 = delta_r0;
    delta_phi1 = delta_phi0;
    z1 = z0;
  }
  if (!valid0 && !valid1) { return false;
}

  double delta_r = 0.0;
  double delta_phi = 0.0;
  if (valid0 && valid1 && !same_sample)
  {
    delta_r = delta_r0 + phi_fraction * (delta_r1 - delta_r0);
    const double unwrapped_delta_phi1 = unwrap_phi_near(delta_phi1, delta_phi0);
    delta_phi = delta_phi0 + phi_fraction * (unwrapped_delta_phi1 - delta_phi0);
    z = z0 + phi_fraction * (z1 - z0);
  }
  else if (valid0)
  {
    delta_r = delta_r0;
    delta_phi = delta_phi0;
    z = z0;
  }
  else
  {
    delta_r = delta_r1;
    delta_phi = delta_phi1;
    z = z1;
  }

  const double radius = hit_radius + delta_r;
  const double output_phi = unwrapped_hit_phi + delta_phi;
  x = radius * std::cos(output_phi);
  y = radius * std::sin(output_phi);
  return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

bool TpcCrossingFinder::make_xyz_point(TrkrDefs::hitsetkey hsk,
                                       TrkrDefs::hitkey hk,
                                       short crossing,
                                       Point& p) const
{
  const unsigned int layer = TrkrDefs::getLayer(hsk);
  const unsigned int hit_side = TpcDefs::getSide(hsk);
  const unsigned int pad = TpcDefs::getPad(hk);
  const unsigned int tbin = TpcDefs::getTBin(hk);

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  if (!sample_drift_lookup(layer, hit_side, pad, tbin, crossing, x, y, z)) { return false;
}

  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.layer = layer;
  p.side = hit_side;
  p.pad = pad;
  p.tbin = tbin;
  p.x = x;
  p.y = y;
  p.z = z;
  return true;
}

bool TpcCrossingFinder::find_time_extrema(const Tpc_AssembledTrack* track,
                                          TrkrDefs::hitsetkey& min_hsk,
                                          TrkrDefs::hitkey& min_hk,
                                          TrkrDefs::hitsetkey& max_hsk,
                                          TrkrDefs::hitkey& max_hk) const
{
  if (!track || track->size_hit_indices() == 0) { return false;
}

  unsigned int min_tbin = std::numeric_limits<unsigned int>::max();
  unsigned int max_tbin = 0;
  bool found = false;
  for (unsigned int ih = 0; ih < track->size_hit_indices(); ++ih)
  {
    const Tpc_AssembledTrack::HitIndex hi = track->get_hit_index(ih);
    const unsigned int tbin = TpcDefs::getTBin(hi.second);
    if (!found || tbin < min_tbin)
    {
      min_tbin = tbin;
      min_hsk = hi.first;
      min_hk = hi.second;
    }
    if (!found || tbin > max_tbin)
    {
      max_tbin = tbin;
      max_hsk = hi.first;
      max_hk = hi.second;
    }
    found = true;
  }
  return found;
}

std::set<short> TpcCrossingFinder::get_available_crossings() const
{
  std::set<short> available_crossings;
  if (!m_vertexMap) { return available_crossings;
}

  for (auto & iter : *m_vertexMap)
  {
    const SvtxVertex* vertex = iter.second;
    if (!vertex) { continue;
}
    available_crossings.insert(vertex->get_beam_crossing());
  }
  return available_crossings;
}

std::set<short> TpcCrossingFinder::get_intt_crossings() const
{
  std::set<short> intt_crossings;
  if (!m_clusterMap) { return intt_crossings;
}

  for (const auto& hitsetkey : m_clusterMap->getHitSetKeys())
  {
    if (static_cast<TrkrDefs::TrkrId>(TrkrDefs::getTrkrId(hitsetkey)) != TrkrDefs::TrkrId::inttId) { continue;
}
    intt_crossings.insert(static_cast<short>(InttDefs::getTimeBucketId(hitsetkey)));
  }
  return intt_crossings;
}

std::map<short, std::vector<TpcCrossingFinder::SiliconVertexHypothesis>>
TpcCrossingFinder::get_vertices_by_crossing() const
{
  std::map<short, std::vector<SiliconVertexHypothesis>> vertices_by_crossing;
  if (!m_vertexMap) { return vertices_by_crossing;
}

  for (auto & iter : *m_vertexMap)
  {
    const SvtxVertex* vertex = iter.second;
    if (!vertex) { continue;
}

    SiliconVertexHypothesis hyp;
    hyp.crossing = vertex->get_beam_crossing();
    hyp.vertex_id = iter.first;
    hyp.x = vertex->get_x();
    hyp.y = vertex->get_y();
    hyp.z = vertex->get_z();
    const double errzz = vertex->get_error(2, 2);
    hyp.sigma_z = errzz >= 0.0 ? std::sqrt(errzz) : std::numeric_limits<double>::quiet_NaN();
    hyp.ntracks = vertex->size_tracks();
    vertices_by_crossing[hyp.crossing].push_back(hyp);
  }
  return vertices_by_crossing;
}

std::vector<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>>
TpcCrossingFinder::select_representatives(const Tpc_AssembledTrack* track,
                                          TrkrDefs::hitsetkey min_hsk,
                                          TrkrDefs::hitkey min_hk,
                                          TrkrDefs::hitsetkey max_hsk,
                                          TrkrDefs::hitkey max_hk) const
{
  std::vector<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>> representatives;
  auto add_unique = [&representatives](const Tpc_AssembledTrack::HitIndex& hi)
  {
    if (std::find(representatives.begin(), representatives.end(), hi) == representatives.end())
    {
      representatives.push_back(hi);
    }
  };

  add_unique({min_hsk, min_hk});
  add_unique({max_hsk, max_hk});
  if (!track) { return representatives;
}

  std::map<unsigned int, Tpc_AssembledTrack::HitIndex> hit_by_layer;
  for (unsigned int ih = 0; ih < track->size_hit_indices(); ++ih)
  {
    const Tpc_AssembledTrack::HitIndex hi = track->get_hit_index(ih);
    hit_by_layer.emplace(TrkrDefs::getLayer(hi.first), hi);
  }
  if (hit_by_layer.empty()) { return representatives;
}

  add_unique(hit_by_layer.begin()->second);
  add_unique(hit_by_layer.rbegin()->second);

  std::vector<unsigned int> layers;
  layers.reserve(hit_by_layer.size());
  for (const auto& item : hit_by_layer) { layers.push_back(item.first);
}
  const unsigned int fractions[] = {1U, 2U};
  for (const unsigned int fraction : fractions)
  {
    if (layers.size() < 3U) { continue;
}
    const std::size_t index = (fraction * (layers.size() - 1U)) / 3U;
    add_unique(hit_by_layer[layers[index]]);
  }

  return representatives;
}

TpcCrossingFinder::ZFitResult TpcCrossingFinder::estimate_tpc_z0_diagnostics(std::vector<Point> points) const
{
  ZFitResult result;
  if (points.size() < 2U) { return result;
}
  std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
    return a.layer < b.layer;
  });

  std::vector<double> s(points.size(), 0.0);
  for (std::size_t i = 1; i < points.size(); ++i)
  {
    s[i] = s[i - 1] + std::hypot(points[i].x - points[i - 1].x, points[i].y - points[i - 1].y);
  }

  std::vector<Tpc_FittingTools::FitPoint> fit_points;
  fit_points.reserve(points.size());
  std::vector<Tpc_FittingTools::FitPoint> radial_fit_points;
  radial_fit_points.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    fit_points.emplace_back(s[i], points[i].z, 1.0);
    radial_fit_points.emplace_back(std::hypot(points[i].x, points[i].y), points[i].z, 1.0);
  }
  const Tpc_FittingTools::LineFit zfit = Tpc_FittingTools::fitLine(fit_points);
  const Tpc_FittingTools::LineFit radial_zfit = Tpc_FittingTools::fitLine(radial_fit_points);
  if (!zfit.ok || !radial_zfit.ok) { return result;
}

  double best_s = s.front();
  double best_d2 = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i + 1U < points.size(); ++i)
  {
    const double vx = points[i + 1U].x - points[i].x;
    const double vy = points[i + 1U].y - points[i].y;
    const double len2 = vx * vx + vy * vy;
    double frac = 0.0;
    if (len2 > 0.0)
    {
      frac = clamp_unit(-(points[i].x * vx + points[i].y * vy) / len2);
    }
    const double x = points[i].x + frac * vx;
    const double y = points[i].y + frac * vy;
    const double d2 = x * x + y * y;
    if (d2 < best_d2)
    {
      best_d2 = d2;
      best_s = s[i] + frac * (s[i + 1U] - s[i]);
    }
  }

  result.valid = std::isfinite(best_s) && std::isfinite(best_d2);
  result.slope = zfit.slope;
  result.intercept = zfit.intercept;
  result.chi2 = zfit.chi2;
  result.ndf = zfit.ndof;
  result.s_at_pca = best_s;
  result.minimum_radius = std::sqrt(best_d2);
  result.z_at_pca = zfit.slope * best_s + zfit.intercept;
  result.z_at_r0 = radial_zfit.intercept;
  result.points = std::move(points);
  result.path_length.reserve(s.size());
  for (const double value : s) { result.path_length.push_back(finite_float_or_nan(value));
}
  result.valid = result.valid && std::isfinite(result.z_at_pca) && std::isfinite(result.z_at_r0);
  return result;
}

bool TpcCrossingFinder::estimate_tpc_z0(std::vector<Point>& points, double& z0) const
{
  const ZFitResult result = estimate_tpc_z0_diagnostics(points);
  if (!result.valid) { return false;
}
  z0 = result.z_at_r0;
  return true;
}

bool TpcCrossingFinder::point_in_tpc(const Point& p) const
{
  const double r = std::hypot(p.x, p.y);
  const double inner = m_idealPadMap ? m_idealPadMap->get_radius(FirstLayer) : 0.0;
  const double outer = m_idealPadMap ? m_idealPadMap->get_radius(LastLayer) : 0.0;
  return std::isfinite(r) && std::isfinite(p.z) &&
    r >= inner - m_radialTolerance &&
    r <= outer + m_radialTolerance &&
    p.z >= -m_tpcHalfLength - m_zTolerance &&
    p.z <= m_tpcHalfLength + m_zTolerance;
}

bool TpcCrossingFinder::point_in_correct_side(const Point& p) const
{
  if (p.side == 0U) { return p.z >= -m_tpcHalfLength - m_zTolerance && p.z <= m_centralMembraneTolerance;
}
  if (p.side == 1U) { return p.z >= -m_centralMembraneTolerance && p.z <= m_tpcHalfLength + m_zTolerance;
}
  return false;
}

TpcCrossingFinder::Candidate
TpcCrossingFinder::test_candidate(const Tpc_AssembledTrack* track,
                                  short crossing,
                                  TrkrDefs::hitsetkey min_hsk,
                                  TrkrDefs::hitkey min_hk,
                                  TrkrDefs::hitsetkey max_hsk,
                                  TrkrDefs::hitkey max_hk,
                                  const std::map<short, std::vector<SiliconVertexHypothesis>>& vertices_by_crossing) const
{
  Candidate candidate;
  candidate.crossing = crossing;
  candidate.silicon_vertex_id = std::numeric_limits<unsigned int>::max();
  candidate.tpc_z0 = std::numeric_limits<double>::quiet_NaN();
  candidate.silicon_vertex_z = std::numeric_limits<double>::quiet_NaN();
  candidate.delta_z = std::numeric_limits<double>::quiet_NaN();

  TpcCrossingCandidate& qa = candidate.qa;
  qa.crossing = crossing;
  qa.is_available_from_intt = true;
  qa.passes_time_window = true;
  qa.was_tested = true;
  qa.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::Unknown);
  qa.max_lookup_time_ns = finite_float_or_nan(m_maxLookupTimeNs);
  qa.min_tbin_time_crossing0_ns = finite_float_or_nan((static_cast<double>(TpcDefs::getTBin(min_hk)) - m_t0) * m_tpcAdcClock);
  qa.max_tbin_time_crossing0_ns = finite_float_or_nan((static_cast<double>(TpcDefs::getTBin(max_hk)) - m_t0) * m_tpcAdcClock);
  qa.candidate_min_time_ns = finite_float_or_nan(static_cast<double>(qa.min_tbin_time_crossing0_ns) - static_cast<double>(crossing) * m_crossingPeriodNs);
  qa.candidate_max_time_ns = finite_float_or_nan(static_cast<double>(qa.max_tbin_time_crossing0_ns) - static_cast<double>(crossing) * m_crossingPeriodNs);
  qa.min_time_margin_ns = qa.candidate_min_time_ns;
  qa.max_time_margin_ns = finite_float_or_nan(m_maxLookupTimeNs - static_cast<double>(qa.candidate_max_time_ns));
  qa.candidate_qa_bits |= FromInttCrossing | PassedTimeWindow;

  auto fill_point_summary = [](const Point& point,
                               unsigned int& layer,
                               unsigned int& side,
                               unsigned int& pad,
                               unsigned int& tbin,
                               float& x,
                               float& y,
                               float& z,
                               float& r,
                               float& phi)
  {
    layer = point.layer;
    side = point.side;
    pad = point.pad;
    tbin = point.tbin;
    x = finite_float_or_nan(point.x);
    y = finite_float_or_nan(point.y);
    z = finite_float_or_nan(point.z);
    r = finite_float_or_nan(std::hypot(point.x, point.y));
    phi = finite_float_or_nan(std::atan2(point.y, point.x));
  };

  auto distance_to_padplane = [this](const Point& point) -> float
  {
    const double padplane_z = point.side == 0U ? -m_tpcHalfLength : m_tpcHalfLength;
    return finite_float_or_nan(std::abs(point.z - padplane_z));
  };

  Point min_point;
  Point max_point;
  qa.min_endpoint_garfield_valid = make_xyz_point(min_hsk, min_hk, crossing, min_point);
  qa.max_endpoint_garfield_valid = make_xyz_point(max_hsk, max_hk, crossing, max_point);
  if (qa.min_endpoint_garfield_valid)
  {
    fill_point_summary(min_point, qa.min_time_layer, qa.min_time_side, qa.min_time_pad, qa.min_time_tbin,
                       qa.min_time_x, qa.min_time_y, qa.min_time_z, qa.min_time_r, qa.min_time_phi);
    qa.min_time_distance_to_central_membrane = finite_float_or_nan(std::abs(min_point.z));
    qa.min_time_distance_to_padplane = distance_to_padplane(min_point);
    qa.min_time_inside_tpc = point_in_tpc(min_point);
    qa.min_time_correct_side = point_in_correct_side(min_point);
    qa.candidate_qa_bits |= MinEndpointGarfieldOK;
  }
  if (qa.max_endpoint_garfield_valid)
  {
    fill_point_summary(max_point, qa.max_time_layer, qa.max_time_side, qa.max_time_pad, qa.max_time_tbin,
                       qa.max_time_x, qa.max_time_y, qa.max_time_z, qa.max_time_r, qa.max_time_phi);
    qa.max_time_distance_to_central_membrane = finite_float_or_nan(std::abs(max_point.z));
    qa.max_time_distance_to_padplane = distance_to_padplane(max_point);
    qa.max_time_inside_tpc = point_in_tpc(max_point);
    qa.max_time_correct_side = point_in_correct_side(max_point);
    qa.candidate_qa_bits |= MaxEndpointGarfieldOK;
  }

  if (!qa.min_endpoint_garfield_valid || !qa.max_endpoint_garfield_valid)
  {
    qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::GarfieldValid);
    candidate.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::GarfieldInvalid);
    qa.rejection_status = candidate.rejection_status;
    return candidate;
  }

  qa.endpoints_inside_tpc = qa.min_time_inside_tpc && qa.max_time_inside_tpc;
  if (qa.endpoints_inside_tpc) { qa.candidate_qa_bits |= EndpointsInsideTPC;
}
  if (!qa.endpoints_inside_tpc)
  {
    qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::InsideTpc);
    candidate.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::OutsideTpcVolume);
    qa.rejection_status = candidate.rejection_status;
    return candidate;
  }

  qa.endpoints_on_correct_side = qa.min_time_correct_side && qa.max_time_correct_side;
  if (qa.endpoints_on_correct_side) { qa.candidate_qa_bits |= EndpointsCorrectSide;
}
  if (!qa.endpoints_on_correct_side)
  {
    qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::CorrectSide);
    candidate.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::WrongTpcSide);
    qa.rejection_status = candidate.rejection_status;
    return candidate;
  }

  std::vector<Point> points;
  for (const auto& hi : select_representatives(track, min_hsk, min_hk, max_hsk, max_hk))
  {
    Point p;
    if (!make_xyz_point(hi.first, hi.second, crossing, p)) { continue;
}
    if (!point_in_tpc(p) || !point_in_correct_side(p)) { continue;
}
    points.push_back(p);
  }

  const ZFitResult fit = estimate_tpc_z0_diagnostics(points);
  qa.fit_valid = fit.valid;
  if (!fit.valid)
  {
    qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::FitValid);
    candidate.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::FitFailed);
    qa.rejection_status = candidate.rejection_status;
    return candidate;
  }

  qa.candidate_qa_bits |= FitSuccessful;
  qa.n_fit_points = fit.points.size();
  qa.z_vs_s_slope = finite_float_or_nan(fit.slope);
  qa.z_vs_s_intercept = finite_float_or_nan(fit.intercept);
  qa.z_fit_chi2 = finite_float_or_nan(fit.chi2);
  qa.z_fit_ndf = fit.ndf;
  qa.s_at_transverse_pca = finite_float_or_nan(fit.s_at_pca);
  qa.minimum_transverse_radius = finite_float_or_nan(fit.minimum_radius);
  qa.tpc_z_at_transverse_pca = finite_float_or_nan(fit.z_at_pca);
  qa.tpc_z_at_r0 = finite_float_or_nan(fit.z_at_r0);
  for (std::size_t i = 0; i < fit.points.size(); ++i)
  {
    const Point& point = fit.points[i];
    qa.fit_point_layer.push_back(point.layer);
    qa.fit_point_pad.push_back(point.pad);
    qa.fit_point_tbin.push_back(point.tbin);
    qa.fit_point_x.push_back(finite_float_or_nan(point.x));
    qa.fit_point_y.push_back(finite_float_or_nan(point.y));
    qa.fit_point_z.push_back(finite_float_or_nan(point.z));
    qa.fit_point_r.push_back(finite_float_or_nan(std::hypot(point.x, point.y)));
    qa.fit_point_s.push_back(i < fit.path_length.size() ? fit.path_length[i] : TpcCrossingCandidate::nan());
  }

  std::map<unsigned int, Tpc_AssembledTrack::HitIndex> hit_by_layer;
  if (track)
  {
    for (unsigned int ih = 0; ih < track->size_hit_indices(); ++ih)
    {
      const Tpc_AssembledTrack::HitIndex hi = track->get_hit_index(ih);
      hit_by_layer.emplace(TrkrDefs::getLayer(hi.first), hi);
    }
  }
  if (!hit_by_layer.empty())
  {
    Point inner;
    if (make_xyz_point(hit_by_layer.begin()->second.first, hit_by_layer.begin()->second.second, crossing, inner))
    {
      qa.inner_layer = inner.layer;
      qa.inner_tbin = inner.tbin;
      qa.inner_x = finite_float_or_nan(inner.x);
      qa.inner_y = finite_float_or_nan(inner.y);
      qa.inner_z = finite_float_or_nan(inner.z);
    }
    Point outer;
    if (make_xyz_point(hit_by_layer.rbegin()->second.first, hit_by_layer.rbegin()->second.second, crossing, outer))
    {
      qa.outer_layer = outer.layer;
      qa.outer_tbin = outer.tbin;
      qa.outer_x = finite_float_or_nan(outer.x);
      qa.outer_y = finite_float_or_nan(outer.y);
      qa.outer_z = finite_float_or_nan(outer.z);
    }
  }

  candidate.tpc_z0 = fit.z_at_r0;
  if (!std::isfinite(candidate.tpc_z0) || std::abs(candidate.tpc_z0 - m_collisionZ) > m_maxCandidateVertexZ)
  {
    qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::VertexCompatible);
    candidate.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::VertexIncompatible);
    qa.rejection_status = candidate.rejection_status;
    return candidate;
  }

  const auto vertex_iter = vertices_by_crossing.find(crossing);
  if (vertex_iter != vertices_by_crossing.end())
  {
    qa.n_vertices_at_crossing = vertex_iter->second.size();
    double best_abs_delta = std::numeric_limits<double>::max();
    for (const SiliconVertexHypothesis& vertex : vertex_iter->second)
    {
      const double delta_z = fit.z_at_r0 - vertex.z;
      const double abs_delta = std::abs(delta_z);
      const double pull_z = std::isfinite(vertex.sigma_z) && vertex.sigma_z > 0.0 ? delta_z / vertex.sigma_z : std::numeric_limits<double>::quiet_NaN();
      qa.candidate_vertex_ids.push_back(vertex.vertex_id);
      qa.candidate_vertex_z.push_back(finite_float_or_nan(vertex.z));
      qa.candidate_vertex_sigma_z.push_back(finite_float_or_nan(vertex.sigma_z));
      qa.candidate_vertex_ntracks.push_back(vertex.ntracks);
      qa.candidate_vertex_delta_z.push_back(finite_float_or_nan(delta_z));
      qa.candidate_vertex_pull_z.push_back(finite_float_or_nan(pull_z));
      if (!std::isfinite(abs_delta) || abs_delta >= best_abs_delta) { continue;
}

      best_abs_delta = abs_delta;
      candidate.has_silicon_vertex = true;
      candidate.silicon_vertex_id = vertex.vertex_id;
      candidate.silicon_vertex_z = vertex.z;
      candidate.delta_z = delta_z;
      qa.closest_vertex_id = vertex.vertex_id;
      qa.closest_vertex_z = finite_float_or_nan(vertex.z);
      qa.closest_vertex_sigma_z = finite_float_or_nan(vertex.sigma_z);
      qa.closest_vertex_delta_z = finite_float_or_nan(delta_z);
      qa.closest_vertex_abs_delta_z = finite_float_or_nan(abs_delta);
      qa.closest_vertex_pull_z = finite_float_or_nan(pull_z);
    }
  }

  qa.has_silicon_vertex = candidate.has_silicon_vertex;
  if (qa.has_silicon_vertex) { qa.candidate_qa_bits |= HasSiliconVertex;
}
  candidate.vertex_compatible = candidate.has_silicon_vertex && std::abs(candidate.delta_z) <= m_maxVertexDz;
  qa.vertex_compatible = candidate.vertex_compatible;
  if (qa.vertex_compatible) { qa.candidate_qa_bits |= PassesVertexDz;
}

  candidate.tpc_valid = true;
  qa.tpc_valid = candidate.tpc_valid;
  qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::FinalRanking);
  qa.rejection_status = candidate.rejection_status;
  return candidate;
}

int TpcCrossingFinder::process_event(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) { return Fun4AllReturnCodes::ABORTEVENT;
}
  if (!m_assembledTracks || !m_decisions) { return Fun4AllReturnCodes::EVENT_OK;
}
  m_decisions->Reset();

  std::set<short> available_crossings = get_available_crossings();
  const std::set<short> intt_crossings = get_intt_crossings();
  available_crossings.insert(intt_crossings.begin(), intt_crossings.end());
  const auto vertices_by_crossing = get_vertices_by_crossing();
  const bool has_vertex_map = m_vertexMap != nullptr;

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event << " available candidate crossings:";
    for (const short crossing : available_crossings) { std::cout << " " << crossing;
}
    std::cout << " | intt crossings:";
    for (const short crossing : intt_crossings) { std::cout << " " << crossing;
}
    std::cout << " | vertices by crossing:";
    for (const auto& item : vertices_by_crossing)
    {
      std::cout << " " << item.first << "(" << item.second.size() << ":";
      for (std::size_t ivtx = 0; ivtx < item.second.size(); ++ivtx)
      {
        const SiliconVertexHypothesis& vertex = item.second[ivtx];
        if (ivtx != 0U) { std::cout << ",";
}
        std::cout << " id=" << vertex.vertex_id
                  << " x=" << vertex.x
                  << " y=" << vertex.y
                  << " z=" << vertex.z;
      }
      std::cout << ")";
    }
    std::cout << std::endl;
  }

  const unsigned int nassembled = m_assembledTracks->size();
  for (unsigned int iassembled = 0; iassembled < nassembled; ++iassembled)
  {
    const Tpc_AssembledTrack* assembled = m_assembledTracks->get_track(iassembled);
    if (!assembled) { continue;
}

    TpcCrossingDecisionv1* decision = new TpcCrossingDecisionv1();
    std::vector<TpcCrossingCandidate> candidate_qa_records;
    auto add_decision_with_candidates = [&candidate_qa_records, decision, this]()
    {
      for (const TpcCrossingCandidate& candidate : candidate_qa_records) { decision->add_candidate(candidate);
}
      m_decisions->add_decision(decision);
    };
    decision->set_assembled_track_id(assembled->get_track_id());
    decision->set_number_of_available_crossings(static_cast<unsigned short>(std::min<std::size_t>(available_crossings.size(), std::numeric_limits<unsigned short>::max())));

    if (available_crossings.empty())
    {
      decision->set_status(TpcCrossingStatus::NoInttCrossings);
      add_decision_with_candidates();
      continue;
    }

    TrkrDefs::hitsetkey min_hsk = 0;
    TrkrDefs::hitkey min_hk = 0;
    TrkrDefs::hitsetkey max_hsk = 0;
    TrkrDefs::hitkey max_hk = 0;
    if (!find_time_extrema(assembled, min_hsk, min_hk, max_hsk, max_hk))
    {
      for (const short crossing : available_crossings)
      {
        TpcCrossingCandidate qa;
        qa.crossing = crossing;
        qa.is_available_from_intt = intt_crossings.contains(crossing);
        qa.candidate_qa_bits = qa.is_available_from_intt ? FromInttCrossing : 0U;
        qa.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::NoValidCrossing);
        qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::PassedTimeWindow);
        candidate_qa_records.push_back(qa);
      }
      decision->set_status(TpcCrossingStatus::NoValidCrossing);
      add_decision_with_candidates();
      continue;
    }

    const unsigned int min_tbin = TpcDefs::getTBin(min_hk);
    const unsigned int max_tbin = TpcDefs::getTBin(max_hk);
    std::vector<short> allowed_crossings;
    for (const short crossing : available_crossings)
    {
      const double min_time0_ns = (static_cast<double>(min_tbin) - m_t0) * m_tpcAdcClock;
      const double max_time0_ns = (static_cast<double>(max_tbin) - m_t0) * m_tpcAdcClock;
      const double min_time_ns = min_time0_ns - static_cast<double>(crossing) * m_crossingPeriodNs;
      const double max_time_ns = max_time0_ns - static_cast<double>(crossing) * m_crossingPeriodNs;
      const bool passes_time_window = std::isfinite(min_time_ns) && std::isfinite(max_time_ns) &&
        min_time_ns > 0.0 && max_time_ns <= m_maxLookupTimeNs;
      if (passes_time_window)
      {
        allowed_crossings.push_back(crossing);
      }
      else
      {
        TpcCrossingCandidate qa;
        qa.crossing = crossing;
        qa.is_available_from_intt = intt_crossings.contains(crossing);
        qa.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::NoAllowedCrossing);
        qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::PassedTimeWindow);
        qa.min_tbin_time_crossing0_ns = finite_float_or_nan(min_time0_ns);
        qa.max_tbin_time_crossing0_ns = finite_float_or_nan(max_time0_ns);
        qa.candidate_min_time_ns = finite_float_or_nan(min_time_ns);
        qa.candidate_max_time_ns = finite_float_or_nan(max_time_ns);
        qa.max_lookup_time_ns = finite_float_or_nan(m_maxLookupTimeNs);
        qa.min_time_margin_ns = finite_float_or_nan(min_time_ns);
        qa.max_time_margin_ns = finite_float_or_nan(m_maxLookupTimeNs - max_time_ns);
        qa.min_time_layer = TrkrDefs::getLayer(min_hsk);
        qa.min_time_side = TpcDefs::getSide(min_hsk);
        qa.min_time_pad = TpcDefs::getPad(min_hk);
        qa.min_time_tbin = min_tbin;
        qa.max_time_layer = TrkrDefs::getLayer(max_hsk);
        qa.max_time_side = TpcDefs::getSide(max_hsk);
        qa.max_time_pad = TpcDefs::getPad(max_hk);
        qa.max_time_tbin = max_tbin;
        qa.candidate_qa_bits = qa.is_available_from_intt ? FromInttCrossing : 0U;
        candidate_qa_records.push_back(qa);
      }
    }
    decision->set_number_of_allowed_crossings(static_cast<unsigned short>(std::min<std::size_t>(allowed_crossings.size(), std::numeric_limits<unsigned short>::max())));

    if (allowed_crossings.empty())
    {
      decision->set_status(TpcCrossingStatus::NoAllowedCrossing);
      add_decision_with_candidates();
      continue;
    }

    std::vector<Candidate> valid_candidates;
    unsigned char last_rejection = static_cast<unsigned char>(TpcCrossingStatus::NoValidCrossing);
    for (const short crossing : allowed_crossings)
    {
      Candidate candidate = test_candidate(assembled, crossing, min_hsk, min_hk, max_hsk, max_hk, vertices_by_crossing);
      decision->set_number_of_tested_crossings(decision->get_number_of_tested_crossings() + 1U);
      candidate.qa.is_available_from_intt = intt_crossings.contains(crossing);
      if (!candidate.qa.is_available_from_intt) { candidate.qa.candidate_qa_bits &= ~FromInttCrossing;
}
      candidate_qa_records.push_back(candidate.qa);
      if (candidate.tpc_valid)
      {
        valid_candidates.push_back(candidate);
      }
      else if (candidate.rejection_status != 0U)
      {
        last_rejection = candidate.rejection_status;
      }
    }
    for (Candidate& candidate : valid_candidates)
    {
      if (candidate.has_silicon_vertex)
      {
        const double abs_delta = std::isfinite(candidate.delta_z) ? std::abs(candidate.delta_z) : std::numeric_limits<double>::max();
        candidate.confidence_tier = abs_delta <= m_maxVertexDz ? 0U : 1U;
        candidate.confidence_score = abs_delta;
      }
      else
      {
        const double abs_z0 = std::isfinite(candidate.tpc_z0) ? std::abs(candidate.tpc_z0 - m_collisionZ) : std::numeric_limits<double>::max();
        candidate.confidence_tier = 2U;
        candidate.confidence_score = abs_z0;
        if (abs_z0 > m_maxTier2BeamlineZ)
        {
          candidate.tpc_valid = false;
          candidate.rejection_status = static_cast<unsigned char>(TpcCrossingStatus::VertexIncompatible);
          last_rejection = candidate.rejection_status;
        }
      }
    }

    std::vector<Candidate> ranked_candidates;
    for (const Candidate& candidate : valid_candidates)
    {
      if (!candidate.tpc_valid) { continue;
}
      if (candidate.confidence_tier == InvalidConfidenceTier) { continue;
}
      if (!std::isfinite(candidate.confidence_score)) { continue;
}
      ranked_candidates.push_back(candidate);
    }

    for (TpcCrossingCandidate& qa : candidate_qa_records)
    {
      const auto iq = std::find_if(valid_candidates.begin(), valid_candidates.end(), [&qa](const Candidate& candidate) { return candidate.crossing == qa.crossing; });
      if (iq == valid_candidates.end()) { continue;
}
      qa.confidence_tier = iq->confidence_tier;
      qa.confidence_score = finite_float_or_nan(iq->confidence_score);
      qa.tpc_valid = iq->tpc_valid;
      if (!iq->tpc_valid && iq->rejection_status != 0U) { qa.rejection_status = iq->rejection_status;
}
    }

    decision->set_number_of_tpc_valid_crossings(static_cast<unsigned short>(std::min<std::size_t>(ranked_candidates.size(), std::numeric_limits<unsigned short>::max())));
    const auto n_vertex_compatible = std::count_if(ranked_candidates.begin(), ranked_candidates.end(), [](const Candidate& candidate) { return candidate.confidence_tier == 0U; });
    decision->set_number_of_vertex_compatible_crossings(static_cast<unsigned short>(std::min<std::size_t>(n_vertex_compatible, std::numeric_limits<unsigned short>::max())));

    if (ranked_candidates.empty())
    {
      decision->set_status(last_rejection != 0U ? last_rejection : static_cast<unsigned char>(TpcCrossingStatus::NoValidCrossing));
      add_decision_with_candidates();
      continue;
    }

    std::sort(ranked_candidates.begin(), ranked_candidates.end(), [](const Candidate& a, const Candidate& b) {
      if (a.confidence_tier != b.confidence_tier) { return a.confidence_tier < b.confidence_tier;
}
      if (a.confidence_score != b.confidence_score) { return a.confidence_score < b.confidence_score;
}
      return a.crossing < b.crossing;
    });

    const Candidate& best_candidate = ranked_candidates.front();
    const double best_score = best_candidate.confidence_score;
    double second_score = std::numeric_limits<double>::infinity();
    for (std::size_t icandidate = 1; icandidate < ranked_candidates.size(); ++icandidate)
    {
      if (ranked_candidates[icandidate].confidence_tier != best_candidate.confidence_tier) { break;
}
      second_score = ranked_candidates[icandidate].confidence_score;
      break;
    }
    const bool close_rival = std::isfinite(second_score) && second_score - best_score <= m_minBestSecondSeparation;

    auto selected_status = TpcCrossingStatus::SelectedByContainment;
    if (best_candidate.confidence_tier == 0U)
    {
      selected_status = close_rival ? TpcCrossingStatus::SelectedByVertexAmbiguous : TpcCrossingStatus::SelectedByVertex;
    }
    else if (best_candidate.confidence_tier == 1U)
    {
      selected_status = TpcCrossingStatus::SelectedByVertexLoose;
    }
    else if (best_candidate.confidence_tier == 2U)
    {
      selected_status = close_rival ? TpcCrossingStatus::SelectedByContainmentAmbiguous : TpcCrossingStatus::SelectedByContainment;
    }

    decision->set_best_abs_delta_z(finite_float_or_nan(best_score));
    decision->set_second_best_abs_delta_z(finite_float_or_nan(second_score));
    decision->set_selected_crossing(best_candidate.crossing);
    decision->set_tpc_z0(finite_float_or_nan(best_candidate.tpc_z0));
    decision->set_silicon_vertex_id(best_candidate.silicon_vertex_id);
    decision->set_silicon_vertex_z(finite_float_or_nan(best_candidate.silicon_vertex_z));
    decision->set_delta_z(finite_float_or_nan(best_candidate.delta_z));
    decision->set_selected_tier(best_candidate.confidence_tier);
    decision->set_selected_score(finite_float_or_nan(best_candidate.confidence_score));
    decision->set_status(selected_status);

    std::vector<short> by_vertex;
    by_vertex.reserve(valid_candidates.size());
    std::vector<short> by_collision;
    by_collision.reserve(valid_candidates.size());
    for (const Candidate& candidate : valid_candidates)
    {
      by_vertex.push_back(candidate.crossing);
      by_collision.push_back(candidate.crossing);
    }
    std::sort(by_vertex.begin(), by_vertex.end(), [&valid_candidates](short a, short b) {
      const auto ia = std::find_if(valid_candidates.begin(), valid_candidates.end(), [a](const Candidate& c) { return c.crossing == a; });
      const auto ib = std::find_if(valid_candidates.begin(), valid_candidates.end(), [b](const Candidate& c) { return c.crossing == b; });
      const double da = std::isfinite(ia->delta_z) ? std::abs(ia->delta_z) : std::numeric_limits<double>::max();
      const double db = std::isfinite(ib->delta_z) ? std::abs(ib->delta_z) : std::numeric_limits<double>::max();
      if (da != db) { return da < db;
}
      return a < b;
    });
    std::sort(by_collision.begin(), by_collision.end(), [this, &valid_candidates](short a, short b) {
      const auto ia = std::find_if(valid_candidates.begin(), valid_candidates.end(), [a](const Candidate& c) { return c.crossing == a; });
      const auto ib = std::find_if(valid_candidates.begin(), valid_candidates.end(), [b](const Candidate& c) { return c.crossing == b; });
      const double da = std::isfinite(ia->tpc_z0) ? std::abs(ia->tpc_z0 - m_collisionZ) : std::numeric_limits<double>::max();
      const double db = std::isfinite(ib->tpc_z0) ? std::abs(ib->tpc_z0 - m_collisionZ) : std::numeric_limits<double>::max();
      if (da != db) { return da < db;
}
      return a < b;
    });

    for (TpcCrossingCandidate& qa : candidate_qa_records)
    {
      const auto iv = std::find(by_vertex.begin(), by_vertex.end(), qa.crossing);
      if (iv != by_vertex.end()) { qa.candidate_rank_by_abs_delta_z = std::distance(by_vertex.begin(), iv);
}
      const auto ic = std::find(by_collision.begin(), by_collision.end(), qa.crossing);
      if (ic != by_collision.end()) { qa.candidate_rank_by_abs_collision_z = std::distance(by_collision.begin(), ic);
}
      const auto iq = std::find_if(valid_candidates.begin(), valid_candidates.end(), [&qa](const Candidate& candidate) { return candidate.crossing == qa.crossing; });
      if (iq != valid_candidates.end())
      {
        qa.confidence_tier = iq->confidence_tier;
        qa.confidence_score = finite_float_or_nan(iq->confidence_score);
        qa.tpc_valid = iq->tpc_valid;
        if (!iq->tpc_valid && iq->rejection_status != 0U) { qa.rejection_status = iq->rejection_status;
}
      }
      qa.is_selected = qa.crossing == decision->get_selected_crossing();
      if (qa.is_selected)
      {
        qa.candidate_qa_bits |= IsSelected;
        qa.first_failed_stage = static_cast<unsigned char>(TpcCrossingCandidateStage::Selected);
      }
      if (qa.candidate_rank_by_abs_delta_z == 0) { qa.candidate_qa_bits |= IsBestVertexCandidate;
}
    }

    add_decision_with_candidates();
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " assembled_tracks=" << nassembled
              << " decisions=" << m_decisions->size()
              << " available_crossings=" << available_crossings.size()
              << " vertex_map=" << has_vertex_map
              << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
