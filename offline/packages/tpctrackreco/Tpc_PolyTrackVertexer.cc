#include "Tpc_PolyTrackVertexer.h"

#include "Tpc_PolyTrack.h"
#include "Tpc_PolyTrackContainer.h"
#include "Tpc_PolyTrackVertexContainerv1.h"
#include "Tpc_PolyTrackVertexv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    while (phi > TMath::Pi()) phi -= 2.0 * TMath::Pi();
    while (phi <= -TMath::Pi()) phi += 2.0 * TMath::Pi();
    return phi;
  }

  bool good(double x)
  {
    return std::isfinite(x) && std::fabs(x) < 1.0e30;
  }
}

Tpc_PolyTrackVertexer::Tpc_PolyTrackVertexer(const std::string& name)
  : SubsysReco(name)
  , m_inputNodeName("TPC_POLYTRACKS")
  , m_outputNodeName("TPC_POLYTRACKVERTICES")
{
}

int Tpc_PolyTrackVertexer::InitRun(PHCompositeNode* topNode)
{
  if (!createNodes(topNode)) return Fun4AllReturnCodes::ABORTRUN;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool Tpc_PolyTrackVertexer::getNodes(PHCompositeNode* topNode)
{
  m_polyTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_inputNodeName.c_str());
  if (!m_polyTracks)
  {
    const char* candidate_names[] = {
      "TPC_POLYTRACKS",
      "Tpc_PolyTracks",
      "Tpc_PolyTrackContainer",
      "TPC_POLYTRACKCONTAINER"
    };

    for (unsigned int i = 0;
         i < sizeof(candidate_names) / sizeof(candidate_names[0]) && !m_polyTracks;
         ++i)
    {
      m_polyTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, candidate_names[i]);
      if (m_polyTracks) m_inputNodeName = candidate_names[i];
    }
  }

  if (!m_polyTracks)
  {
    std::cerr << Name() << "::getNodes - missing Tpc_PolyTrackContainer node" << std::endl;
    return false;
  }

  return true;
}

bool Tpc_PolyTrackVertexer::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_vertices = findNode::getClass<Tpc_PolyTrackVertexContainer>(topNode, m_outputNodeName.c_str());
  if (!m_vertices)
  {
    m_vertices = new Tpc_PolyTrackVertexContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_vertices, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return true;
}

Tpc_PolyTrackVertexer::TrackVertexFit
Tpc_PolyTrackVertexer::fitTrack(const Tpc_PolyTrack* trk) const
{
  TrackVertexFit fit;
  if (!trk || trk->get_fit_status() <= 0) return fit;

  const double x = trk->get_x();
  const double y = trk->get_y();
  const double z = trk->get_z();
  const double px = trk->get_px();
  const double py = trk->get_py();
  const double pz = trk->get_pz();
  const double charge = trk->get_charge();
  if (!good(x) || !good(y) || !good(z) ||
      !good(px) || !good(py) || !good(pz) || !good(charge))
  {
    return fit;
  }

  const double pt2 = px * px + py * py;
  if (pt2 <= 1.0e-20) return fit;
  const double pt = std::sqrt(pt2);

  if (std::fabs(charge * m_magneticFieldTesla) < 1.0e-12)
  {
    const double scale = -(x * px + y * py) / pt2;
    fit.pca_x = x + scale * px;
    fit.pca_y = y + scale * py;
    fit.pca_z = z + scale * pz;
    fit.dca2d = (x * py - y * px) / pt;
  }
  else
  {
    const double signed_radius = pt / (0.003 * charge * m_magneticFieldTesla);
    const double radius = std::fabs(signed_radius);
    if (!good(radius) || radius <= 0.0) return fit;

    const double tx = px / pt;
    const double ty = py / pt;
    const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
    const double xc = x + sign * radius * ty;
    const double yc = y - sign * radius * tx;
    const double dc = std::hypot(xc, yc);
    if (!good(dc) || dc <= 1.0e-12) return fit;

    fit.pca_x = xc * (1.0 - radius / dc);
    fit.pca_y = yc * (1.0 - radius / dc);

    const double phi0 = std::atan2(y - yc, x - xc);
    const double pca_phi_on_circle = std::atan2(fit.pca_y - yc, fit.pca_x - xc);
    double best_arc = 0.0;
    double best_abs_arc = 1.0e30;
    for (int k = -4; k <= 4; ++k)
    {
      const double dphi = pca_phi_on_circle - phi0 + 2.0 * TMath::Pi() * static_cast<double>(k);
      const double arc = -sign * radius * dphi;
      const double abs_arc = std::fabs(arc);
      if (abs_arc < best_abs_arc)
      {
        best_abs_arc = abs_arc;
        best_arc = arc;
      }
    }

    fit.pca_z = z + (pz / pt) * best_arc;
    fit.dca2d = dc - radius;
  }

  fit.pca_radius = std::hypot(fit.pca_x, fit.pca_y);
  fit.pca_phi = wrap_phi(std::atan2(fit.pca_y, fit.pca_x));
  fit.z0 = fit.pca_z;
  fit.track_id = trk->get_track_id();
  fit.source_assembled_track_id = trk->get_source_assembled_track_id();
  fit.nclusters = trk->get_nclusters();
  fit.ok = good(fit.dca2d) && good(fit.z0) && good(fit.pca_radius) && good(fit.pca_phi);
  return fit;
}

Tpc_PolyTrackVertexer::CollisionFit
Tpc_PolyTrackVertexer::fitCollision(const std::vector<TrackVertexFit>& tracks) const
{
  CollisionFit fit;
  if (tracks.size() < 2) return fit;

  double sum_w = 0.0;
  double sum_x = 0.0;
  double sum_y = 0.0;
  double sum_z = 0.0;
  for (const TrackVertexFit& trk : tracks)
  {
    if (!trk.ok) continue;
    const double w = std::max(1.0, static_cast<double>(trk.nclusters));
    sum_w += w;
    sum_x += w * trk.pca_x;
    sum_y += w * trk.pca_y;
    sum_z += w * trk.z0;
  }

  if (sum_w <= 0.0) return fit;
  fit.x = sum_x / sum_w;
  fit.y = sum_y / sum_w;
  fit.z = sum_z / sum_w;

  double sum_res2 = 0.0;
  unsigned int nvalid = 0;
  for (const TrackVertexFit& trk : tracks)
  {
    if (!trk.ok) continue;
    const double w = std::max(1.0, static_cast<double>(trk.nclusters));
    const double dz = trk.z0 - fit.z;
    sum_res2 += w * dz * dz;
    ++nvalid;
  }

  if (nvalid < 2) return fit;
  fit.z_rms = std::sqrt(sum_res2 / sum_w);
  fit.ntracks = nvalid;
  fit.ok = good(fit.x) && good(fit.y) && good(fit.z) && good(fit.z_rms);
  return fit;
}

std::vector<Tpc_PolyTrackVertexer::CollisionFit>
Tpc_PolyTrackVertexer::fitCollisions(std::vector<TrackVertexFit> tracks) const
{
  std::vector<CollisionFit> collisions;
  tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                              [](const TrackVertexFit& trk)
                              {
                                return !trk.ok;
                              }),
               tracks.end());
  if (tracks.size() < 2) return collisions;

  std::sort(tracks.begin(), tracks.end(),
            [](const TrackVertexFit& a, const TrackVertexFit& b)
            {
              return a.z0 < b.z0;
            });

  const double separation = std::max(0.0, m_collisionZSeparation);
  std::vector<TrackVertexFit> cluster;
  cluster.reserve(tracks.size());
  double cluster_sum_w = 0.0;
  double cluster_sum_z = 0.0;

  for (const TrackVertexFit& trk : tracks)
  {
    const double w = std::max(1.0, static_cast<double>(trk.nclusters));
    const double cluster_z = cluster_sum_w > 0.0 ? cluster_sum_z / cluster_sum_w : trk.z0;

    if (!cluster.empty() && trk.z0 - cluster_z >= separation)
    {
      const CollisionFit collision = fitCollision(cluster);
      if (collision.ok) collisions.push_back(collision);
      cluster.clear();
      cluster_sum_w = 0.0;
      cluster_sum_z = 0.0;
    }

    cluster.push_back(trk);
    cluster_sum_w += w;
    cluster_sum_z += w * trk.z0;
  }

  const CollisionFit collision = fitCollision(cluster);
  if (collision.ok) collisions.push_back(collision);
  return collisions;
}

int Tpc_PolyTrackVertexer::process_event(PHCompositeNode* topNode)
{
  if (!m_polyTracks || !m_vertices)
  {
    if (!getNodes(topNode) || !createNodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;
  }

  m_vertices->Reset();

  std::vector<TrackVertexFit> collision_tracks;
  const unsigned int ntracks = m_polyTracks ? m_polyTracks->size() : 0;
  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const Tpc_PolyTrack* trk = m_polyTracks->get_track(itrk);
    const TrackVertexFit fit = fitTrack(trk);
    if (!fit.ok) continue;

    Tpc_PolyTrackVertexv1* out = new Tpc_PolyTrackVertexv1();
    out->set_track_id(fit.track_id);
    out->set_source_assembled_track_id(fit.source_assembled_track_id);
    out->set_dca2d(fit.dca2d);
    out->set_z0(fit.z0);
    out->set_pca_valid(1);
    out->set_pca_x(fit.pca_x);
    out->set_pca_y(fit.pca_y);
    out->set_pca_z(fit.pca_z);
    out->set_pca_radius(fit.pca_radius);
    out->set_pca_phi(fit.pca_phi);
    m_vertices->add_vertex(out);

    if (fit.nclusters >= m_collisionMinClusters) collision_tracks.push_back(fit);
  }

  const std::vector<CollisionFit> collisions = fitCollisions(collision_tracks);
  m_vertices->set_collision_min_clusters(m_collisionMinClusters);
  m_vertices->clear_collision_vertices();
  for (const CollisionFit& collision : collisions)
  {
    m_vertices->add_collision_vertex(collision.x, collision.y, collision.z,
                                     collision.z_rms, collision.ntracks);
  }
  m_vertices->set_collision_vertex_valid(collisions.empty() ? 0 : 1);

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - input poly_tracks=" << ntracks
              << " vertices=" << m_vertices->size()
              << " collision_vertices=" << m_vertices->get_collision_vertex_count()
              << " first_collision_ntracks=" << m_vertices->get_collision_ntracks()
              << " first_collision_x=" << m_vertices->get_collision_x()
              << " first_collision_y=" << m_vertices->get_collision_y()
              << " first_collision_z=" << m_vertices->get_collision_z()
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
