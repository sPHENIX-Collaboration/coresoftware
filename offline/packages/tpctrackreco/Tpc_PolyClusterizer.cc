#include "Tpc_PolyClusterizer.h"

#include "IdealPadMap.h"
#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackContainer.h"
#include "Tpc_PolyClusterContainerv1.h"
#include "Tpc_PolyClusterv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TPolyLine3D.h>
#include <phgarfield/PHGarfield.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace
{
  constexpr unsigned int FirstLayer = 7;
  constexpr unsigned int LastLayer = 54;
  constexpr unsigned int NLayers = LastLayer - FirstLayer + 1;
  constexpr unsigned int NSides = 2;
  constexpr unsigned int NSectors = 12;
  constexpr double PhiConsistencyTolerance = 1.0e-10;

  double wrap_phi(double phi)
  {
    while (phi > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi <= -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }

  double unwrap_phi_near(double phi, const double reference)
  {
    while (phi - reference > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi - reference <= -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }

  double clamp_unit(const double value)
  {
    return std::max(0.0, std::min(1.0, value));
  }

  double phi_sample_fraction(const unsigned int sample)
  {
    if (Tpc_PolyClusterizer::NPhiSamples <= 1U)
    {
      return 0.0;
    }
    return static_cast<double>(sample) / static_cast<double>(Tpc_PolyClusterizer::NPhiSamples - 1U);
  }
}  // namespace

Tpc_PolyClusterizer::Tpc_PolyClusterizer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("TPC_ASSEMBLEDTRACKS")
  , m_outputNodeName("TPC_POLYCLUSTERS")
{
}

Tpc_PolyClusterizer::~Tpc_PolyClusterizer()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
  delete m_garfield;
  m_garfield = nullptr;
}

int Tpc_PolyClusterizer::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - failed to load IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  delete m_garfield;
  // m_garfield = new PHGarfield(Name() + "_PHGarfield");

  const std::string electricFieldMap = "/sphenix/user/mitrankov/garf/include/sphenix_rossegger_garfield_field.root";
  // sphenix_3d_ibf_field_new.root sphenix_rossegger_garfield_field.root;

  m_garfield = new PHGarfield(Name() + "_PHGarfield", electricFieldMap, m_kEffSide0, m_kEffSide1);
  configure_garfield(m_garfield);
  if (m_garfield->InitRun(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cerr << Name() << "::InitRun - PHGarfield InitRun failed" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (!build_drift_lookup())
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

void Tpc_PolyClusterizer::configure_garfield(PHGarfield* garfield) const
{
  if (!garfield)
  {
    return;
  }

  garfield->MoveTpc(m_tpcMove[0], m_tpcMove[1], m_tpcMove[2]);
  for (const auto& rotation : m_tpcRotations)
  {
    garfield->RotateTpc(rotation[0], rotation[1], rotation[2]);
  }
  garfield->SetCMVoltageDefault(m_cmVoltageDefault);
}

int Tpc_PolyClusterizer::getNodes(PHCompositeNode* topNode)
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

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyClusterizer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_outputNodeName);
  if (!m_clusters)
  {
    m_clusters = new Tpc_PolyClusterContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_clusters, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

unsigned int Tpc_PolyClusterizer::drift_lookup_index(const unsigned int layer_index,
                                                     const unsigned int side,
                                                     const unsigned int sector,
                                                     const unsigned int sample)
{
  return (((layer_index * NSides + side) * NSectors + sector) * NPhiSamples + sample);
}

bool Tpc_PolyClusterizer::build_drift_lookup()
{
  if (!m_idealPadMap || !m_garfield)
  {
    return false;
  }
  if (m_reverseDriftStepNs <= 0.0 || !std::isfinite(m_reverseDriftStepNs))
  {
    return false;
  }

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
    if (!std::isfinite(radius) || pads_per_sector == 0U)
    {
      std::cerr << Name() << "::build_drift_lookup - invalid geometry for layer " << layer << std::endl;
      return false;
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
          if (!std::isfinite(phi_local) || !std::isfinite(phi_global))
          {
            std::cerr << Name() << "::build_drift_lookup - invalid phi for layer " << layer
                      << " side " << side
                      << " sector " << sector
                      << " local_phibin " << local_phibin
                      << " global_pad " << global_pad
                      << " phi_local " << phi_local
                      << " phi_global " << phi_global << std::endl;
            return false;
          }

          const double phi_difference = wrap_phi(phi_global - phi_local);
          if (std::abs(phi_difference) > PhiConsistencyTolerance)
          {
            std::cerr << Name() << "::build_drift_lookup - inconsistent IdealPadMap phi overloads"
                      << " layer " << layer
                      << " side " << side
                      << " sector " << sector
                      << " local_phibin " << local_phibin
                      << " global_pad " << global_pad
                      << " phi_local " << phi_local
                      << " phi_global " << phi_global
                      << " phi_difference " << phi_difference << std::endl;
            return false;
          }

          const double phi = phi_local;
          const double x0 = radius * std::cos(phi);
          const double y0 = radius * std::sin(phi);
          TPolyLine3D* drift = m_garfield->ReverseDrift(x0, y0, z0, m_reverseDriftStepNs);
          if (!drift || drift->GetN() <= 0)
          {
            delete drift;
            std::cerr << Name() << "::build_drift_lookup - ReverseDrift failed for layer " << layer
                      << " side " << side << " sector " << sector << " sample " << sample << std::endl;
            return false;
          }

          const int npoints = drift->GetN();
          const Float_t* xyz = drift->GetP();
          if (!xyz || npoints <= 0)
          {
            delete drift;
            std::cerr << Name() << "::build_drift_lookup - empty drift points for layer " << layer
                      << " side " << side << " sector " << sector << " sample " << sample << std::endl;
            return false;
          }

          DriftPolyline& polyline = m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
          polyline.phi = phi;
          polyline.points.resize(static_cast<std::size_t>(npoints));
          for (int ipoint = 0; ipoint < npoints; ++ipoint)
          {
            const int idx = 3 * ipoint;
            DriftPoint& point = polyline.points[static_cast<std::size_t>(ipoint)];
            const double output_r = std::hypot(static_cast<double>(xyz[idx]), static_cast<double>(xyz[idx + 1]));
            const double output_phi = unwrap_phi_near(std::atan2(static_cast<double>(xyz[idx + 1]), static_cast<double>(xyz[idx])), phi);
            point.delta_r = static_cast<float>(output_r - radius);
            point.delta_phi = static_cast<float>(output_phi - phi);
            point.z = xyz[idx + 2];
          }
          delete drift;
          ++nbuilt;
        }
      }
    }
  }

  for (unsigned int layer = FirstLayer; layer <= LastLayer; ++layer)
  {
    const unsigned int layer_index = layer - FirstLayer;
    const double radius = m_idealPadMap->get_radius(layer);
    for (unsigned int side = 0; side < NSides; ++side)
    {
      const DriftPolyline& reference = m_driftLookup[drift_lookup_index(layer_index, side, 0, 0)];
      bool symmetry_ok = !reference.points.empty() && std::isfinite(radius);
      unsigned int length_mismatches = 0;
      double max_delta_r = 0.0;
      double max_delta_phi_arc = 0.0;
      double max_delta_z = 0.0;

      for (unsigned int sector = 0; sector < NSectors; ++sector)
      {
        for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
        {
          const DriftPolyline& polyline = m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
          if (polyline.points.size() != reference.points.size())
          {
            ++length_mismatches;
            symmetry_ok = false;
          }

          const std::size_t npoints = std::min(reference.points.size(), polyline.points.size());
          for (std::size_t ipoint = 0; ipoint < npoints; ++ipoint)
          {
            const DriftPoint& ref_point = reference.points[ipoint];
            const DriftPoint& point = polyline.points[ipoint];
            max_delta_r = std::max(max_delta_r, std::abs(static_cast<double>(point.delta_r - ref_point.delta_r)));
            max_delta_phi_arc = std::max(max_delta_phi_arc,
                                         std::abs(radius * wrap_phi(static_cast<double>(point.delta_phi - ref_point.delta_phi))));
            max_delta_z = std::max(max_delta_z, std::abs(static_cast<double>(point.z - ref_point.z)));
          }
        }
      }

      if (max_delta_r > 0.001 ||
          max_delta_phi_arc > 0.001 ||
          max_delta_z > 0.001)
      {
        symmetry_ok = false;
      }

      std::cout << Name() << "::build_drift_lookup - phi symmetry "
                << (symmetry_ok ? "holds" : "broken")
                << " layer=" << layer
                << " side=" << side
                << " max_dr=" << max_delta_r
                << " max_r_dphi=" << max_delta_phi_arc
                << " max_dz=" << max_delta_z
                << " length_mismatches=" << length_mismatches << std::endl;
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::build_drift_lookup - built " << nbuilt << " drift polylines" << std::endl;
  }
  return nbuilt == NLayers * NSides * NSectors * NPhiSamples;
}

bool Tpc_PolyClusterizer::sample_drift_lookup(const unsigned int layer,
                                              const unsigned int side,
                                              const unsigned int pad,
                                              const unsigned int tbin,
                                              double& x,
                                              double& y,
                                              double& z) const
{
  if (!m_idealPadMap)
  {
    return false;
  }
  if (layer < FirstLayer || layer > LastLayer)
  {
    return false;
  }
  if (side >= NSides)
  {
    return false;
  }
  if (m_reverseDriftStepNs <= 0.0 || !std::isfinite(m_reverseDriftStepNs))
  {
    return false;
  }

  const unsigned int pads_per_sector = m_idealPadMap->get_pads_per_sector_for_layer(layer);
  if (pads_per_sector == 0U)
  {
    return false;
  }

  const unsigned int sector = pad / pads_per_sector;
  if (sector >= NSectors)
  {
    return false;
  }

  const double hit_radius = m_idealPadMap->get_radius(layer);
  const double hit_phi = m_idealPadMap->get_phi(side, layer, pad);
  if (!std::isfinite(hit_radius) || !std::isfinite(hit_phi))
  {
    return false;
  }

  const double corrected_tbin = static_cast<double>(tbin) - m_t0;
  const double target_time_ns = corrected_tbin * m_tpcAdcClock;
  if (target_time_ns <= 0.0 || !std::isfinite(target_time_ns))
  {
    return false;
  }

  const unsigned int layer_index = layer - FirstLayer;
  std::array<const DriftPolyline*, NPhiSamples> samples{};
  std::array<double, NPhiSamples> sample_phi{};
  for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
  {
    samples[sample] = &m_driftLookup[drift_lookup_index(layer_index, side, sector, sample)];
    if (samples[sample]->points.empty())
    {
      return false;
    }
    sample_phi[sample] = samples[sample]->phi;
    if (sample > 0U)
    {
      sample_phi[sample] = unwrap_phi_near(sample_phi[sample], sample_phi[sample - 1U]);
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
    sample0 = 0;
    sample1 = 0;
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
      const bool in_interval = increasing ? (unwrapped_hit_phi >= sample_phi[sample] && unwrapped_hit_phi <= sample_phi[sample + 1U]) : (unwrapped_hit_phi <= sample_phi[sample] && unwrapped_hit_phi >= sample_phi[sample + 1U]);
      if (!in_interval)
      {
        continue;
      }

      sample0 = sample;
      sample1 = sample + 1U;
      const double denom = sample_phi[sample1] - sample_phi[sample0];
      phi_fraction = (denom != 0.0) ? clamp_unit((unwrapped_hit_phi - sample_phi[sample0]) / denom) : 0.0;
      bracket_found = true;
      break;
    }
  }

  if (!bracket_found)
  {
    std::cerr << Name() << "::sample_drift_lookup - failed to bracket hit phi"
              << " layer " << layer
              << " side " << side
              << " sector " << sector
              << " pad " << pad
              << " hit_phi " << hit_phi
              << " unwrapped_hit_phi " << unwrapped_hit_phi
              << " sample_phi";
    for (unsigned int sample = 0; sample < NPhiSamples; ++sample)
    {
      std::cerr << " " << sample_phi[sample];
    }
    std::cerr << std::endl;
    return false;
  }

  auto sample_time = [this, target_time_ns](const DriftPolyline& polyline,
                                            double& delta_r,
                                            double& delta_phi,
                                            double& point_z) -> bool
  {
    const int npoints = static_cast<int>(polyline.points.size());
    if (npoints <= 0)
    {
      return false;
    }

    const double max_time_ns = static_cast<double>(npoints - 1) * m_reverseDriftStepNs;
    if (target_time_ns > max_time_ns)
    {
      return false;
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

  if (!valid0 && !valid1)
  {
    return false;
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

bool Tpc_PolyClusterizer::make_xyz_point(TrkrDefs::hitsetkey hsk,
                                         TrkrDefs::hitkey hk,
                                         Point& p) const
{
  if (!m_hits || !m_idealPadMap)
  {
    return false;
  }

  TrkrHitSet* hitset = m_hits->findHitSet(hsk);
  if (!hitset)
  {
    return false;
  }
  TrkrHit* hit = hitset->getHit(hk);
  if (!hit)
  {
    return false;
  }

  const unsigned int layer = TrkrDefs::getLayer(hsk);
  const unsigned int hit_side = TpcDefs::getSide(hsk);
  const unsigned int pad = TpcDefs::getPad(hk);
  const unsigned int tbin = TpcDefs::getTBin(hk);
  const double adc = hit->getAdc();

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  if (!sample_drift_lookup(layer, hit_side, pad, tbin, x, y, z))
  {
    return false;
  }

  p.hitsetkey = hsk;
  p.hitkey = hk;
  p.layer = layer;
  p.side = hit_side;
  p.pad = pad;
  p.tbin = tbin;
  p.adc = adc;
  p.x = x;
  p.y = y;
  p.z = z;
  return true;
}

Tpc_PolyClusterizer::ClusterParameters
Tpc_PolyClusterizer::make_cluster_parameters(const std::vector<Point>& points,
                                             const Centroid& centroid,
                                             const int side) const
{
  ClusterParameters params;
  if (points.empty() || !centroid.ok || !m_idealPadMap)
  {
    return params;
  }

  std::set<unsigned int> pads;
  std::set<unsigned int> tbins;
  std::map<unsigned int, double> adc_by_pad;
  for (const Point& p : points)
  {
    params.adc += p.adc;
    pads.insert(p.pad);
    tbins.insert(p.tbin);
    adc_by_pad[p.pad] += p.adc;
  }

  params.phi_width = static_cast<unsigned int>(pads.size());
  params.time_width = static_cast<unsigned int>(tbins.size());

  unsigned int max_adc_pad = 0;
  double max_adc = -std::numeric_limits<double>::max();
  for (const auto& pad_adc : adc_by_pad)
  {
    if (pad_adc.second > max_adc)
    {
      max_adc = pad_adc.second;
      max_adc_pad = pad_adc.first;
    }
  }

  const unsigned int total_phibins = m_idealPadMap->get_total_phibins(centroid.layer);
  const double pad_phi_width = total_phibins > 0U ? 2.0 * M_PI / static_cast<double>(total_phibins) : 0.0;
  const double cluster_phi = std::atan2(centroid.y, centroid.x);
  const double max_adc_phi = m_idealPadMap->get_phi(static_cast<unsigned int>(side), centroid.layer, max_adc_pad);
  if (pad_phi_width > 0.0 && std::isfinite(cluster_phi) && std::isfinite(max_adc_phi))
  {
    params.phase = wrap_phi(cluster_phi - max_adc_phi) / pad_phi_width;
  }

  return params;
}

Tpc_PolyClusterizer::Centroid
Tpc_PolyClusterizer::make_centroid(const std::vector<Point>& points)
{
  Centroid c;
  if (points.empty())
  {
    return c;
  }

  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  for (const Point& p : points)
  {
    sx += p.x;
    sy += p.y;
    sz += p.z;
  }

  const double n = static_cast<double>(points.size());
  c.x = sx / n;
  c.y = sy / n;
  c.z = sz / n;

  double sxx = 0.0;
  double syy = 0.0;
  double szz = 0.0;
  for (const Point& p : points)
  {
    const double dx = p.x - c.x;
    const double dy = p.y - c.y;
    const double dz = p.z - c.z;
    sxx += dx * dx;
    syy += dy * dy;
    szz += dz * dz;
  }

  c.rms_x = std::sqrt(sxx / n);
  c.rms_y = std::sqrt(syy / n);
  c.rms_z = std::sqrt(szz / n);
  c.layer = points.front().layer;
  c.ok = std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z);
  return c;
}

int Tpc_PolyClusterizer::process_event(PHCompositeNode* /*unused*/)
{
  if (!m_assembledTracks || !m_clusters || !m_garfield)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  m_clusters->Reset();

  const unsigned int nassembled = m_assembledTracks->size();
  unsigned int nclusters = 0;

  for (int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      for (unsigned int iassembled = 0; iassembled < nassembled; ++iassembled)
      {
        const Tpc_AssembledTrack* assembled = m_assembledTracks->get_track(iassembled);
        if (!assembled)
        {
          continue;
        }
        if (assembled->get_side() != side)
        {
          continue;
        }
        if (assembled->get_first_sector() % 12U != sector)
        {
          continue;
        }

        std::map<unsigned int, std::vector<Point>> points_by_layer;
        for (unsigned int ih = 0; ih < assembled->size_hit_indices(); ++ih)
        {
          const Tpc_AssembledTrack::HitIndex hi = assembled->get_hit_index(ih);
          if (TpcDefs::getSide(hi.first) != static_cast<unsigned int>(side))
          {
            continue;
          }

          Point p;
          if (make_xyz_point(hi.first, hi.second, p))
          {
            points_by_layer[p.layer].push_back(p);
          }
        }
        if (points_by_layer.empty())
        {
          continue;
        }

        for (const auto& layer_points : points_by_layer)
        {
          const std::vector<Point>& points = layer_points.second;
          const Centroid centroid = make_centroid(points);
          if (!centroid.ok)
          {
            continue;
          }

          Tpc_PolyClusterv1* out = new Tpc_PolyClusterv1();
          out->set_event(m_event);
          out->set_cluster_id(m_clusters->size());
          out->set_source_assembled_track_id(assembled->get_track_id());
          out->set_side(side);
          out->set_centroid_x(centroid.x);
          out->set_centroid_y(centroid.y);
          out->set_centroid_z(centroid.z);
          out->set_rms_x(centroid.rms_x);
          out->set_rms_y(centroid.rms_y);
          out->set_rms_z(centroid.rms_z);

          const ClusterParameters params = make_cluster_parameters(points, centroid, static_cast<int>(points.front().side));
          out->set_adc(params.adc);
          out->set_phi_width(params.phi_width);
          out->set_time_width(params.time_width);
          out->set_phase(params.phase);
          for (const Point& p : points)
          {
            out->add_hit(p.hitsetkey, p.hitkey, p.x, p.y, p.z);
          }
          if (out->size_hits() == 0)
          {
            delete out;
            continue;
          }
          m_clusters->add_cluster(out);
          ++nclusters;
        }
      }
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " assembled_tracks=" << nassembled
              << " poly_clusters=" << m_clusters->size()
              << " layer_clusters=" << nclusters << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
