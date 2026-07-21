#include "Tpc_PolyTrackReco.h"

#include "Tpc_PolyTrackContainerv1.h"
#include "Tpc_PolyTrackv1.h"
#include "IdealPadMap.h"
#include "Tpc_PolyCluster.h"
#include "Tpc_PolyClusterContainer.h"
#include "Tpc_FittingTools.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <map>
#include <vector>

Tpc_PolyTrackReco::Tpc_PolyTrackReco(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("TPC_POLYCLUSTERS")
  , m_outputNodeName("TPC_POLYTRACKS")
{
}

Tpc_PolyTrackReco::~Tpc_PolyTrackReco()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;
}

int Tpc_PolyTrackReco::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK) return Fun4AllReturnCodes::ABORTRUN;

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - failed to load IdealPadMap" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

double Tpc_PolyTrackReco::calc_dedx(const std::vector<const Tpc_PolyCluster*>& clusters,
                                          const Tpc_FittingTools::FitResult& fit,
                                          const bool fit_ok) const
{
  if (clusters.empty() || !fit_ok || !m_idealPadMap)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double thickness_per_region[4] = {
    m_idealPadMap->get_layer_thickness(7),
    m_idealPadMap->get_layer_thickness(8),
    m_idealPadMap->get_layer_thickness(27),
    m_idealPadMap->get_layer_thickness(50)
  };

  std::vector<double> dedxlist;
  dedxlist.reserve(clusters.size());
  for (const Tpc_PolyCluster* cluster : clusters)
  {
    if (!cluster || cluster->size_hits() == 0) continue;
    const unsigned int layer = TrkrDefs::getLayer(cluster->get_hit_index(0).first);
    double thick = std::numeric_limits<double>::quiet_NaN();
    if (layer < 23U)
    {
      thick = thickness_per_region[layer % 2U == 0U ? 1 : 0];
    }
    else if (layer < 39U)
    {
      thick = thickness_per_region[2];
    }
    else
    {
      thick = thickness_per_region[3];
    }
    if (!std::isfinite(thick) || thick <= 0.0) continue;

    const double x = cluster->get_centroid_x();
    const double y = cluster->get_centroid_y();
    const double r = std::hypot(x, y);
    if (!std::isfinite(r)) continue;

    double adc = cluster->get_adc() / thick;
    if (!fit.is_line && std::isfinite(fit.curvature))
    {
      const double alpha = 0.5 * r * std::fabs(fit.curvature);
      double alphacorr = std::cos(alpha);
      if (alphacorr < 0.0 || alphacorr > 4.0)
      {
        alphacorr = 4.0;
      }
      adc *= alphacorr;
    }

    adc *= std::clamp(std::sin(fit.theta), 0.0, 4.0);
    if (std::isfinite(adc)) dedxlist.push_back(adc);
  }

  if (dedxlist.empty())
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::sort(dedxlist.begin(), dedxlist.end());
  const unsigned int trunc_max = static_cast<unsigned int>(dedxlist.size() * 0.7);
  double sumdedx = 0.0;
  unsigned int ndedx = 0;
  for (unsigned int j = 0; j <= trunc_max && j < dedxlist.size(); ++j)
  {
    sumdedx += dedxlist[j];
    ++ndedx;
  }

  return ndedx > 0U ? sumdedx / static_cast<double>(ndedx) : std::numeric_limits<double>::quiet_NaN();
}

int Tpc_PolyTrackReco::getNodes(PHCompositeNode* topNode)
{
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_inputNodeName);
  if (!m_clusters)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyTrackReco::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_polyTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_outputNodeName);
  if (!m_polyTracks)
  {
    m_polyTracks = new Tpc_PolyTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_polyTracks, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void Tpc_PolyTrackReco::fillTpc_PolyTrack(unsigned int source_assembled_track_id,
                                             const std::vector<const Tpc_PolyCluster*>& clusters,
                                             const Tpc_FittingTools::FitResult& fit,
                                             const bool fit_ok)
{
  Tpc_PolyTrackv1* out = new Tpc_PolyTrackv1();
  out->set_event(m_event);
  out->set_track_id(m_polyTracks->size());
  out->set_source_assembled_track_id(source_assembled_track_id);
  out->set_fit_status(fit_ok ? 1 : 0);
  out->set_nclusters(static_cast<unsigned int>(clusters.size()));
  out->set_dedx(calc_dedx(clusters, fit, fit_ok));

  if (fit_ok)
  {
    if (fit.is_line)
    {
      out->set_x(fit.line_x);
      out->set_y(fit.line_y);
      out->set_z(fit.line_z);
      out->set_px(fit.line_dx);
      out->set_py(fit.line_dy);
      out->set_pz(fit.line_dz);
      out->set_charge(0.0);
    }
    else
    {
      const double sin_phi = std::sin(fit.phi0);
      const double cos_phi = std::cos(fit.phi0);
      out->set_x(-fit.d0 * sin_phi);
      out->set_y(fit.d0 * cos_phi);
      out->set_z(fit.z0);

      const double abs_curvature = std::fabs(fit.curvature);
      const double pt = abs_curvature > 0.0 ? 0.003 * (m_magneticFieldTesla) / abs_curvature : 0.0;
      const double tan_theta = std::tan(fit.theta);
      out->set_px(pt * cos_phi);
      out->set_py(pt * sin_phi);
      const double pz = std::fabs(tan_theta) > 1.0e-12 ? (pt / tan_theta) : 0.0;
      out->set_pz(pz);
      double charge = fit.curvature >= 0.0 ? -1.0 : 1.0;
      out->set_charge(charge);
    }
    out->set_chi2(fit.chi2_xy + fit.chi2_z);
    out->set_ndf(static_cast<double>(fit.ndof_xy + fit.ndof_z));
  }

  m_polyTracks->add_track(out);
}

int Tpc_PolyTrackReco::process_event(PHCompositeNode* topNode)
{
  if (!m_clusters || !m_polyTracks)
  {
    if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK ||
        createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
  }

  m_polyTracks->Reset();

  std::map<unsigned int, std::vector<const Tpc_PolyCluster*>> clusters_by_track;
  const unsigned int nclusters = m_clusters->size();
  for (unsigned int icluster = 0; icluster < nclusters; ++icluster)
  {
    const Tpc_PolyCluster* cluster = m_clusters->get_cluster(icluster);
    if (!cluster) continue;
    clusters_by_track[cluster->get_source_assembled_track_id()].push_back(cluster);
  }

  for (const auto& track_clusters : clusters_by_track)
  {
    const std::vector<const Tpc_PolyCluster*>& clusters = track_clusters.second;
    std::vector<Tpc_FittingTools::Point> fit_points;
    fit_points.reserve(clusters.size());
    for (const Tpc_PolyCluster* cluster : clusters)
    {
      if (!cluster) continue;
      Tpc_FittingTools::Point fp;
      fp.x = cluster->get_centroid_x();
      fp.y = cluster->get_centroid_y();
      fp.z = cluster->get_centroid_z();
      if (std::isfinite(fp.x) && std::isfinite(fp.y) && std::isfinite(fp.z)) fit_points.push_back(fp);
    }

    Tpc_FittingTools::FitResult fit;
    const bool fit_ok = (m_fitMode == FitMode::Line3D) ?
      Tpc_FittingTools::fitLine3D(fit_points, fit) :
      Tpc_FittingTools::fit(fit_points, fit);
    fillTpc_PolyTrack(track_clusters.first, clusters, fit, fit_ok);
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " poly_clusters=" << nclusters
              << " poly_tracks=" << m_polyTracks->size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}
