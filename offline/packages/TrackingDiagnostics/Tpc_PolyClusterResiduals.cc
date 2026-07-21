#include "Tpc_PolyClusterResiduals.h"

#include "tpctrackreco/Tpc_PolyTrack.h"
#include "tpctrackreco/Tpc_PolyTrackContainer.h"
#include "tpctrackreco/Tpc_PolyTrackVertex.h"
#include "tpctrackreco/Tpc_PolyTrackVertexContainer.h"
#include "tpctrackreco/Tpc_PolyCluster.h"
#include "tpctrackreco/Tpc_PolyClusterContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace
{
  double wrap_phi(double phi)
  {
    const double pi = std::acos(-1.0);
    while (phi > pi) phi -= 2.0 * pi;
    while (phi <= -pi) phi += 2.0 * pi;
    return phi;
  }

  struct HelixCircle
  {
    bool ok {false};
    double xc {0.0};
    double yc {0.0};
    double radius {0.0};
    double sign {0.0};
    double phi0 {0.0};
    double z0 {0.0};
    double dzds {0.0};
  };

  HelixCircle make_track_circle(const Tpc_PolyTrack* trk,
                                const double magnetic_field_tesla)
  {
    HelixCircle circle;
    if (!trk || trk->get_fit_status() == 0) return circle;

    const double x = trk->get_x();
    const double y = trk->get_y();
    const double z = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    const double charge = trk->get_charge();
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(charge))
    {
      return circle;
    }

    const double pt = std::hypot(px, py);
    if (pt <= 0.0 || std::fabs(charge * magnetic_field_tesla) < 1.0e-12) return circle;

    const double signed_radius = pt / (0.003 * charge * magnetic_field_tesla);
    circle.radius = std::fabs(signed_radius);
    if (circle.radius <= 0.0 || !std::isfinite(circle.radius)) return circle;

    const double tx = px / pt;
    const double ty = py / pt;
    circle.sign = signed_radius > 0.0 ? 1.0 : -1.0;
    circle.xc = x + circle.sign * circle.radius * ty;
    circle.yc = y - circle.sign * circle.radius * tx;
    circle.phi0 = std::atan2(y - circle.yc, x - circle.xc);
    circle.z0 = z;
    circle.dzds = pz / pt;
    circle.ok = std::isfinite(circle.xc) && std::isfinite(circle.yc) &&
                std::isfinite(circle.phi0) && std::isfinite(circle.z0) &&
                std::isfinite(circle.dzds);
    return circle;
  }

  bool helix_z_at_dca_to_vertex(const HelixCircle& circle,
                                const double vertex_x,
                                const double vertex_y,
                                double& z_at_dca)
  {
    if (!circle.ok) return false;

    const double dx = vertex_x - circle.xc;
    const double dy = vertex_y - circle.yc;
    const double dc = std::hypot(dx, dy);
    if (!std::isfinite(dc) || dc <= 1.0e-12) return false;

    const double pca_x = circle.xc + circle.radius * dx / dc;
    const double pca_y = circle.yc + circle.radius * dy / dc;
    const double pca_phi_on_circle = std::atan2(pca_y - circle.yc, pca_x - circle.xc);

    double best_arc = 0.0;
    double best_abs_arc = std::numeric_limits<double>::max();
    const double pi = std::acos(-1.0);
    for (int k = -4; k <= 4; ++k)
    {
      const double dphi = pca_phi_on_circle - circle.phi0 + 2.0 * pi * static_cast<double>(k);
      const double arc = -circle.sign * circle.radius * dphi;
      const double abs_arc = std::fabs(arc);
      if (abs_arc < best_abs_arc)
      {
        best_abs_arc = abs_arc;
        best_arc = arc;
      }
    }

    z_at_dca = circle.z0 + circle.dzds * best_arc;
    return std::isfinite(z_at_dca);
  }

  bool helix_z_at_radius(const HelixCircle& circle,
                         const double target_r,
                         const double reference_z,
                         const double arc_direction,
                         double& z_state)
  {
    if (!circle.ok) return false;

    const double center_r = std::hypot(circle.xc, circle.yc);
    if (!std::isfinite(target_r) || !std::isfinite(center_r) ||
        target_r <= 0.0 || center_r <= 1.0e-12)
    {
      return false;
    }

    const double radius_sum = target_r + circle.radius;
    const double radius_diff = std::fabs(target_r - circle.radius);
    if (center_r > radius_sum || center_r < radius_diff) return false;

    const double a = (target_r * target_r - circle.radius * circle.radius + center_r * center_r) /
                     (2.0 * center_r);
    const double h2 = target_r * target_r - a * a;
    if (h2 < -1.0e-8) return false;

    const double h = std::sqrt(std::max(0.0, h2));
    const double ux = circle.xc / center_r;
    const double uy = circle.yc / center_r;

    double best_z = 0.0;
    double best_abs_dz = std::numeric_limits<double>::max();
    const double pi = std::acos(-1.0);
    for (int isign = -1; isign <= 1; isign += 2)
    {
      const double x = a * ux - static_cast<double>(isign) * h * uy;
      const double y = a * uy + static_cast<double>(isign) * h * ux;
      const double phi_on_circle = std::atan2(y - circle.yc, x - circle.xc);

      for (int k = -4; k <= 4; ++k)
      {
        const double dphi = phi_on_circle - circle.phi0 + 2.0 * pi * static_cast<double>(k);
        const double arc = -circle.sign * circle.radius * dphi;
        const double z = circle.z0 + arc_direction * circle.dzds * arc;
        const double abs_dz = std::fabs(z - reference_z);
        if (abs_dz < best_abs_dz)
        {
          best_abs_dz = abs_dz;
          best_z = z;
        }
      }
    }

    if (best_abs_dz == std::numeric_limits<double>::max()) return false;

    z_state = best_z;
    return std::isfinite(z_state);
  }

  bool choose_collision_vertex(const Tpc_PolyTrackVertexContainer* vertices,
                               const HelixCircle& circle,
                               double& vertex_x,
                               double& vertex_y,
                               double& vertex_z)
  {
    if (!vertices || !vertices->get_collision_vertex_valid() || !circle.ok) return false;

    double best_dz = std::numeric_limits<double>::max();
    const unsigned int nvertices = vertices->get_collision_vertex_count();
    for (unsigned int ivtx = 0; ivtx < nvertices; ++ivtx)
    {
      const double x = vertices->get_collision_x(ivtx);
      const double y = vertices->get_collision_y(ivtx);
      const double z = vertices->get_collision_z(ivtx);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;

      double z_at_dca = 0.0;
      if (!helix_z_at_dca_to_vertex(circle, x, y, z_at_dca)) continue;

      const double dz = std::fabs(z - z_at_dca);
      if (dz < best_dz)
      {
        best_dz = dz;
        vertex_x = x;
        vertex_y = y;
        vertex_z = z;
      }
    }

    return best_dz != std::numeric_limits<double>::max();
  }

  bool line_xy_at_z(const Tpc_PolyTrack* trk,
                    const double z,
                    const double arc_direction,
                    double& x_state,
                    double& y_state)
  {
    if (!trk || trk->get_fit_status() == 0 || !std::isfinite(z)) return false;

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        std::fabs(pz) < 1.0e-12)
    {
      return false;
    }

    const double dz = z - z0;
    x_state = x0 + arc_direction * px / pz * dz;
    y_state = y0 + arc_direction * py / pz * dz;
    return std::isfinite(x_state) && std::isfinite(y_state);
  }

  bool line_z_at_radius(const Tpc_PolyTrack* trk,
                        const double target_r,
                        const double reference_z,
                        const double arc_direction,
                        double& z_state)
  {
    if (!trk || trk->get_fit_status() == 0 || !std::isfinite(target_r) || target_r <= 0.0) return false;

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        std::fabs(pz) < 1.0e-12)
    {
      return false;
    }

    const double ax = arc_direction * px / pz;
    const double ay = arc_direction * py / pz;
    const double a = ax * ax + ay * ay;
    const double b = 2.0 * (x0 * ax + y0 * ay);
    const double c = x0 * x0 + y0 * y0 - target_r * target_r;
    if (a < 1.0e-20) return false;

    const double disc = b * b - 4.0 * a * c;
    if (disc < -1.0e-8) return false;
    const double root = std::sqrt(std::max(0.0, disc));
    const double dz1 = (-b - root) / (2.0 * a);
    const double dz2 = (-b + root) / (2.0 * a);
    const double z1 = z0 + dz1;
    const double z2 = z0 + dz2;
    z_state = std::fabs(z1 - reference_z) <= std::fabs(z2 - reference_z) ? z1 : z2;
    return std::isfinite(z_state);
  }

  bool line_z_at_dca_to_vertex(const Tpc_PolyTrack* trk,
                               const double vertex_x,
                               const double vertex_y,
                               const double arc_direction,
                               double& z_at_dca,
                               double& dca_xy)
  {
    if (!trk || trk->get_fit_status() == 0) return false;

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(vertex_x) || !std::isfinite(vertex_y) || std::fabs(pz) < 1.0e-12)
    {
      return false;
    }

    const double ax = arc_direction * px / pz;
    const double ay = arc_direction * py / pz;
    const double den = ax * ax + ay * ay;
    if (den < 1.0e-20) return false;
    const double dz = ((vertex_x - x0) * ax + (vertex_y - y0) * ay) / den;
    const double x_at_dca = x0 + ax * dz;
    const double y_at_dca = y0 + ay * dz;
    z_at_dca = z0 + dz;
    dca_xy = std::hypot(x_at_dca - vertex_x, y_at_dca - vertex_y);
    return std::isfinite(z_at_dca) && std::isfinite(dca_xy);
  }

  bool choose_collision_vertex_line(const Tpc_PolyTrackVertexContainer* vertices,
                                    const Tpc_PolyTrack* trk,
                                    const double arc_direction,
                                    double& vertex_x,
                                    double& vertex_y,
                                    double& vertex_z,
                                    double& rdca)
  {
    if (!vertices || !vertices->get_collision_vertex_valid() || !trk) return false;

    double best_dz = std::numeric_limits<double>::max();
    const unsigned int nvertices = vertices->get_collision_vertex_count();
    for (unsigned int ivtx = 0; ivtx < nvertices; ++ivtx)
    {
      const double x = vertices->get_collision_x(ivtx);
      const double y = vertices->get_collision_y(ivtx);
      const double z = vertices->get_collision_z(ivtx);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;

      double z_at_dca = 0.0;
      double dca_xy = 0.0;
      if (!line_z_at_dca_to_vertex(trk, x, y, arc_direction, z_at_dca, dca_xy)) continue;

      const double dz = std::fabs(z - z_at_dca);
      if (dz < best_dz)
      {
        best_dz = dz;
        vertex_x = x;
        vertex_y = y;
        vertex_z = z;
        rdca = dca_xy;
      }
    }

    return best_dz != std::numeric_limits<double>::max();
  }

  unsigned int cluster_sector(const Tpc_PolyCluster* cluster)
  {
    if (!cluster || cluster->size_hits() == 0) return 0xffffffffu;
    const Tpc_PolyCluster::HitIndex hit_index = cluster->get_hit_index(0);
    return TpcDefs::getSectorId(hit_index.first);
  }

  unsigned int cluster_layer(const Tpc_PolyCluster* cluster)
  {
    if (!cluster || cluster->size_hits() == 0) return 0xffffffffu;
    const Tpc_PolyCluster::HitIndex hit_index = cluster->get_hit_index(0);
    return TrkrDefs::getLayer(hit_index.first);
  }

  bool project_track_to_z(const Tpc_PolyTrack* trk,
                          const double z,
                          const double magnetic_field_tesla,
                          const double arc_direction,
                          const bool use_straight_line,
                          double& x_state,
                          double& y_state)
  {
    if (!trk || trk->get_fit_status() == 0) return false;

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    const double charge = trk->get_charge();

    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(z))
    {
      return false;
    }

    if (use_straight_line || std::fabs(charge * magnetic_field_tesla) < 1.0e-12)
    {
      return line_xy_at_z(trk, z, arc_direction, x_state, y_state);
    }

    const double pt = std::hypot(px, py);
    if (pt <= 0.0 || std::fabs(pz) < 1.0e-12) return false;

    const double signed_radius = pt / (0.003 * charge * magnetic_field_tesla);
    const double radius = std::fabs(signed_radius);
    if (radius <= 0.0 || !std::isfinite(radius)) return false;

    const double tx = px / pt;
    const double ty = py / pt;
    const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
    const double xc = x0 + sign * radius * ty;
    const double yc = y0 - sign * radius * tx;
    const double phi0 = std::atan2(y0 - yc, x0 - xc);
    const double dzds = pz / pt;
    if (std::fabs(dzds) < 1.0e-12) return false;

    const double arc = arc_direction * (z - z0) / dzds;
    const double phi = phi0 - sign * arc / radius;
    x_state = xc + radius * std::cos(phi);
    y_state = yc + radius * std::sin(phi);
    return std::isfinite(x_state) && std::isfinite(y_state);
  }

  double cluster_line_residual2(const Tpc_PolyTrack* poly_track,
                                const std::vector<const Tpc_PolyCluster*>& clusters,
                                const double magnetic_field_tesla,
                                const double arc_direction,
                                const bool use_straight_line)
  {
    if (!poly_track || clusters.empty()) return std::numeric_limits<double>::max();

    double sum = 0.0;
    unsigned int n = 0;
    for (const Tpc_PolyCluster* cluster : clusters)
    {
      if (!cluster || !cluster->isValid()) continue;
      const double cx = cluster->get_centroid_x();
      const double cy = cluster->get_centroid_y();
      const double cz = cluster->get_centroid_z();
      if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cz)) continue;

      double x = 0.0;
      double y = 0.0;
      if (!project_track_to_z(poly_track, cz, magnetic_field_tesla, arc_direction, use_straight_line, x, y)) continue;

      const double dx = x - cx;
      const double dy = y - cy;
      sum += dx * dx + dy * dy;
      ++n;
    }

    return n > 0 ? sum / static_cast<double>(n) : std::numeric_limits<double>::max();
  }

}

Tpc_PolyClusterResiduals::Tpc_PolyClusterResiduals(const std::string& name,
                                                 const std::string& outfilename)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_clusterNodeName("TPC_POLYCLUSTERS")
  , m_finalTrackNodeName("TPC_POLYTRACKS")
  , m_finalTrackVertexNodeName("TPC_POLYTRACKVERTICES")
{
}

Tpc_PolyClusterResiduals::~Tpc_PolyClusterResiduals()
{
  if (m_outfile)
  {
    delete m_outfile;
    m_outfile = nullptr;
  }
}

int Tpc_PolyClusterResiduals::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tree = new TTree("residuals", "TPC poly cluster r-phi residuals");
  m_tree->Branch("event", &m_event, "event/i");
  m_tree->Branch("poly_track_id", &m_finalTrackId, "poly_track_id/i");
  m_tree->Branch("source_cluster_id", &m_sourceClusterId, "source_cluster_id/i");
  m_tree->Branch("source_assembled_track_id", &m_sourceAssembledTrackId, "source_assembled_track_id/i");
  m_tree->Branch("cluster_index", &m_clusterIndex);
  m_tree->Branch("side", &m_side, "side/I");
  m_tree->Branch("sector", &m_sector);
  m_tree->Branch("layer", &m_layer);
  m_tree->Branch("ntpc_clusters", &m_ntpcClusters, "ntpc_clusters/i");
  m_tree->Branch("fit_status", &m_fitStatus, "fit_status/I");
  m_tree->Branch("pt", &m_pt, "pt/D");
  m_tree->Branch("px", &m_px, "px/D");
  m_tree->Branch("py", &m_py, "py/D");
  m_tree->Branch("pz", &m_pz, "pz/D");
  m_tree->Branch("eta", &m_eta, "eta/D");
  m_tree->Branch("theta", &m_theta, "theta/D");
  m_tree->Branch("charge", &m_charge, "charge/D");
  m_tree->Branch("chi2", &m_chi2, "chi2/D");
  m_tree->Branch("ndf", &m_ndf, "ndf/D");
  m_tree->Branch("quality", &m_quality, "quality/D");
  m_tree->Branch("dedx", &m_dedx, "dedx/D");
  m_tree->Branch("vertex_x", &m_vertexX, "vertex_x/D");
  m_tree->Branch("vertex_y", &m_vertexY, "vertex_y/D");
  m_tree->Branch("vertex_z", &m_vertexZ, "vertex_z/D");
  m_tree->Branch("vertex_r", &m_vertexR, "vertex_r/D");
  m_tree->Branch("pca_x", &m_pcaX, "pca_x/D");
  m_tree->Branch("pca_y", &m_pcaY, "pca_y/D");
  m_tree->Branch("pca_z", &m_pcaZ, "pca_z/D");
  m_tree->Branch("rDCA", &m_rDCA, "rDCA/D");
  m_tree->Branch("rDCA_zero", &m_rDCAZero, "rDCA_zero/D");
  m_tree->Branch("zDCA", &m_zDCA, "zDCA/D");
  m_tree->Branch("R", &m_R, "R/D");
  m_tree->Branch("rzslope", &m_rzSlope, "rzslope/D");
  m_tree->Branch("cluster_x", &m_clusterX);
  m_tree->Branch("cluster_y", &m_clusterY);
  m_tree->Branch("cluster_z", &m_clusterZ);
  m_tree->Branch("cluster_r", &m_clusterR);
  m_tree->Branch("cluster_phi", &m_clusterPhi);
  m_tree->Branch("cluster_adc", &m_clusterAdc);
  m_tree->Branch("cluster_pad_size", &m_clusterPadSize);
  m_tree->Branch("state_x", &m_stateX);
  m_tree->Branch("state_y", &m_stateY);
  m_tree->Branch("state_z", &m_stateZ);
  m_tree->Branch("state_z_dca", &m_stateZDca);
  m_tree->Branch("state_r", &m_stateR);
  m_tree->Branch("state_phi", &m_statePhi);
  m_tree->Branch("delta_phi", &m_deltaPhi);
  m_tree->Branch("residual_rphi", &m_residualRPhi);
  m_tree->Branch("residual_z", &m_residualZ);

  return Fun4AllReturnCodes::EVENT_OK;
}

bool Tpc_PolyClusterResiduals::get_nodes(PHCompositeNode* topNode)
{
  m_clusters = findNode::getClass<Tpc_PolyClusterContainer>(topNode, m_clusterNodeName.c_str());
  if (!m_clusters)
  {
    std::cerr << Name() << " - missing " << m_clusterNodeName << std::endl;
    return false;
  }

  m_finalTracks = findNode::getClass<Tpc_PolyTrackContainer>(topNode, m_finalTrackNodeName.c_str());
  if (!m_finalTracks)
  {
    std::cerr << Name() << " - missing " << m_finalTrackNodeName << std::endl;
    return false;
  }

  m_finalTrackVertices = findNode::getClass<Tpc_PolyTrackVertexContainer>(topNode, m_finalTrackVertexNodeName.c_str());
  if (!m_finalTrackVertices && Verbosity() > 0)
  {
    std::cerr << Name() << " - missing " << m_finalTrackVertexNodeName
              << ", rdca will be NaN" << std::endl;
  }

  return true;
}

void Tpc_PolyClusterResiduals::reset_tree_values()
{
  m_event = m_evt;
  m_finalTrackId = 0;
  m_sourceClusterId = 0;
  m_sourceAssembledTrackId = 0;
  m_side = 0;
  m_ntpcClusters = 0;
  m_fitStatus = 0;
  m_pt = 0.0;
  m_px = std::numeric_limits<double>::quiet_NaN();
  m_py = std::numeric_limits<double>::quiet_NaN();
  m_pz = std::numeric_limits<double>::quiet_NaN();
  m_eta = std::numeric_limits<double>::quiet_NaN();
  m_theta = std::numeric_limits<double>::quiet_NaN();
  m_charge = 0.0;
  m_chi2 = 0.0;
  m_ndf = 0.0;
  m_quality = std::numeric_limits<double>::quiet_NaN();
  m_dedx = std::numeric_limits<double>::quiet_NaN();
  m_vertexX = std::numeric_limits<double>::quiet_NaN();
  m_vertexY = std::numeric_limits<double>::quiet_NaN();
  m_vertexZ = std::numeric_limits<double>::quiet_NaN();
  m_vertexR = std::numeric_limits<double>::quiet_NaN();
  m_pcaX = std::numeric_limits<double>::quiet_NaN();
  m_pcaY = std::numeric_limits<double>::quiet_NaN();
  m_pcaZ = std::numeric_limits<double>::quiet_NaN();
  m_zDCA = std::numeric_limits<double>::quiet_NaN();
  m_rDCA = std::numeric_limits<double>::quiet_NaN();
  m_rDCAZero = std::numeric_limits<double>::quiet_NaN();
  m_R = std::numeric_limits<double>::quiet_NaN();
  m_rzSlope = std::numeric_limits<double>::quiet_NaN();
  m_clusterIndex.clear();
  m_sector.clear();
  m_layer.clear();
  m_clusterX.clear();
  m_clusterY.clear();
  m_clusterZ.clear();
  m_clusterR.clear();
  m_clusterPhi.clear();
  m_clusterAdc.clear();
  m_clusterPadSize.clear();
  m_stateX.clear();
  m_stateY.clear();
  m_stateZ.clear();
  m_stateZDca.clear();
  m_stateR.clear();
  m_statePhi.clear();
  m_deltaPhi.clear();
  m_residualRPhi.clear();
  m_residualZ.clear();
}

int Tpc_PolyClusterResiduals::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!get_nodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;
  if (!m_tree) return Fun4AllReturnCodes::EVENT_OK;

  std::map<unsigned int, std::vector<const Tpc_PolyCluster*> > clusters_by_source_assembled_track_id;
  std::map<unsigned int, const Tpc_PolyTrackVertex*> track_vertices_by_track_id;
  std::map<unsigned int, const Tpc_PolyTrackVertex*> track_vertices_by_source_assembled_track_id;
  for (unsigned int icluster = 0; icluster < m_clusters->size(); ++icluster)
  {
    const Tpc_PolyCluster* cluster = m_clusters->get_cluster(icluster);
    if (!cluster || !cluster->isValid()) continue;
    clusters_by_source_assembled_track_id[cluster->get_source_assembled_track_id()].push_back(cluster);
  }

  if (m_finalTrackVertices)
  {
    for (unsigned int ivtx = 0; ivtx < m_finalTrackVertices->size(); ++ivtx)
    {
      const Tpc_PolyTrackVertex* vtx = m_finalTrackVertices->get_vertex(ivtx);
      if (!vtx) continue;
      track_vertices_by_track_id[vtx->get_track_id()] = vtx;
      track_vertices_by_source_assembled_track_id[vtx->get_source_assembled_track_id()] = vtx;
    }
  }

  unsigned int nfilled = 0;
  const unsigned int npoly_tracks = m_finalTracks->size();
  for (unsigned int ifinal = 0; ifinal < npoly_tracks; ++ifinal)
  {
    const Tpc_PolyTrack* poly_track = m_finalTracks->get_track(ifinal);
    if (!poly_track || !poly_track->isValid()) continue;

    const double px = poly_track->get_px();
    const double py = poly_track->get_py();
    const double pz = poly_track->get_pz();
    const double charge = poly_track->get_charge();
    const bool use_straight_line = m_useStraightLineTracks || std::fabs(charge * m_magneticFieldTesla) < 1.0e-12;
    const double pt = std::hypot(px, py);
    if (!std::isfinite(pt) || (!use_straight_line && (pt < m_minPt || pt > m_maxPt))) continue;

    const double eta = (pt > 0.0 && std::isfinite(pz)) ? std::asinh(pz / pt) : std::numeric_limits<double>::quiet_NaN();
    const double theta = (pt > 0.0 && std::isfinite(pz)) ? std::atan2(pt, pz) : std::numeric_limits<double>::quiet_NaN();
    const double chi2 = poly_track->get_chi2();
    const double ndf = poly_track->get_ndf();
    const double quality = (std::isfinite(chi2) && std::isfinite(ndf) && ndf > 0.0) ? chi2 / ndf : std::numeric_limits<double>::quiet_NaN();

    const auto cluster_iter = clusters_by_source_assembled_track_id.find(poly_track->get_source_assembled_track_id());
    if (cluster_iter == clusters_by_source_assembled_track_id.end()) continue;

    const std::vector<const Tpc_PolyCluster*>& track_clusters = cluster_iter->second;
    const Tpc_PolyTrackVertex* track_vertex = nullptr;
    const auto vertex_source_iter = track_vertices_by_source_assembled_track_id.find(poly_track->get_source_assembled_track_id());
    if (vertex_source_iter != track_vertices_by_source_assembled_track_id.end())
    {
      track_vertex = vertex_source_iter->second;
    }
    else
    {
      const auto vertex_track_iter = track_vertices_by_track_id.find(poly_track->get_track_id());
      if (vertex_track_iter != track_vertices_by_track_id.end())
      {
        track_vertex = vertex_track_iter->second;
      }
    }
    const unsigned int ntpc_clusters = track_clusters.size();
    if (ntpc_clusters < m_minTpcClusters || ntpc_clusters > m_maxTpcClusters) continue;

    const double forward_residual2 = cluster_line_residual2(poly_track, track_clusters, m_magneticFieldTesla, 1.0, use_straight_line);
    const double reverse_residual2 = cluster_line_residual2(poly_track, track_clusters, m_magneticFieldTesla, -1.0, use_straight_line);
    const double arc_direction = forward_residual2 <= reverse_residual2 ? 1.0 : -1.0;

    double vertex_x = std::numeric_limits<double>::quiet_NaN();
    double vertex_y = std::numeric_limits<double>::quiet_NaN();
    double vertex_z = std::numeric_limits<double>::quiet_NaN();
    double rdca = std::numeric_limits<double>::quiet_NaN();
    double rdca_zero = std::numeric_limits<double>::quiet_NaN();
    double pca_x = std::numeric_limits<double>::quiet_NaN();
    double pca_y = std::numeric_limits<double>::quiet_NaN();
    double pca_z = std::numeric_limits<double>::quiet_NaN();
    double zdca = std::numeric_limits<double>::quiet_NaN();
    const HelixCircle circle = use_straight_line ? HelixCircle() : make_track_circle(poly_track, m_magneticFieldTesla);
    const double dedx = poly_track->get_dedx();

    if (use_straight_line)
    {
      if (choose_collision_vertex_line(m_finalTrackVertices, poly_track, arc_direction, vertex_x, vertex_y, vertex_z, rdca))
      {
        double z_at_zero = 0.0;
        line_z_at_dca_to_vertex(poly_track, 0.0, 0.0, arc_direction, z_at_zero, rdca_zero);
      }
    }
    else if (choose_collision_vertex(m_finalTrackVertices, circle, vertex_x, vertex_y, vertex_z))
    {
      rdca = std::hypot(circle.xc - vertex_x, circle.yc - vertex_y) - circle.radius;
      rdca_zero = std::hypot(circle.xc, circle.yc) - circle.radius;
    }

    if (track_vertex && track_vertex->get_pca_valid())
    {
      pca_x = track_vertex->get_pca_x();
      pca_y = track_vertex->get_pca_y();
      pca_z = track_vertex->get_pca_z();
      if (std::isfinite(pca_z) && std::isfinite(vertex_z))
      {
        zdca = pca_z - vertex_z;
      }
    }

    reset_tree_values();
    m_event = m_evt;
    m_finalTrackId = poly_track->get_track_id();
    m_sourceClusterId = track_clusters.empty() ? 0 : track_clusters.front()->get_cluster_id();
    m_sourceAssembledTrackId = poly_track->get_source_assembled_track_id();
    m_side = track_clusters.empty() ? 0 : track_clusters.front()->get_side();
    m_ntpcClusters = ntpc_clusters;
    m_fitStatus = poly_track->get_fit_status();
    m_pt = pt;
    m_px = px;
    m_py = py;
    m_pz = pz;
    m_eta = eta;
    m_theta = theta;
    m_charge = charge;
    m_chi2 = chi2;
    m_ndf = ndf;
    m_quality = quality;
    m_dedx = dedx;
    m_vertexX = vertex_x;
    m_vertexY = vertex_y;
    m_vertexZ = vertex_z;
    m_vertexR = std::hypot(vertex_x, vertex_y);
    m_pcaX = pca_x;
    m_pcaY = pca_y;
    m_pcaZ = pca_z;
    m_rDCA = rdca;
    m_rDCAZero = rdca_zero;
    m_zDCA = zdca;
    m_R = circle.ok ? circle.radius : std::numeric_limits<double>::quiet_NaN();
    m_rzSlope = circle.ok ? circle.dzds : ((use_straight_line && pt > 0.0) ? pz / pt : std::numeric_limits<double>::quiet_NaN());
    for (const Tpc_PolyCluster* cluster : track_clusters)
    {
      if (!cluster || !cluster->isValid()) continue;
      const double cluster_x = cluster->get_centroid_x();
      const double cluster_y = cluster->get_centroid_y();
      const double cluster_z = cluster->get_centroid_z();
      if (!std::isfinite(cluster_x) || !std::isfinite(cluster_y) || !std::isfinite(cluster_z)) continue;

      double state_x = 0.0;
      double state_y = 0.0;
      if (!project_track_to_z(poly_track, cluster_z, m_magneticFieldTesla, arc_direction, use_straight_line, state_x, state_y)) continue;

      const int cluster_side = cluster->get_side();
      double state_z = std::numeric_limits<double>::quiet_NaN();
      double residual_z = std::numeric_limits<double>::quiet_NaN();
      const double cluster_r_for_state = std::hypot(cluster_x, cluster_y);
      const bool have_state_z = use_straight_line ?
        line_z_at_radius(poly_track, cluster_r_for_state, cluster_z, arc_direction, state_z) :
        helix_z_at_radius(circle, cluster_r_for_state, cluster_z, arc_direction, state_z);
      if (have_state_z)
      {
        residual_z = cluster_z - state_z;
      }

      m_clusterIndex.push_back(cluster->get_cluster_id());
      m_side = cluster_side;
      m_sector.push_back(cluster_sector(cluster));
      m_layer.push_back(cluster_layer(cluster));
      m_clusterX.push_back(cluster_x);
      m_clusterY.push_back(cluster_y);
      m_clusterZ.push_back(cluster_z);
      const double cluster_r = std::hypot(cluster_x, cluster_y);
      const double cluster_phi = std::atan2(cluster_y, cluster_x);
      m_clusterR.push_back(cluster_r);
      m_clusterPhi.push_back(cluster_phi);
      m_clusterAdc.push_back(cluster->get_adc());
      m_clusterPadSize.push_back(cluster->get_phi_width());
      m_stateX.push_back(state_x);
      m_stateY.push_back(state_y);
      m_stateZ.push_back(state_z);
      m_stateZDca.push_back(state_z);
      m_stateR.push_back(std::hypot(state_x, state_y));
      const double state_phi = std::atan2(state_y, state_x);
      m_statePhi.push_back(state_phi);
      const double delta_phi = wrap_phi(cluster_phi - state_phi);
      m_deltaPhi.push_back(delta_phi);
      m_residualRPhi.push_back(cluster_r * delta_phi);
      m_residualZ.push_back(residual_z);
      ++nfilled;
    }

    if (!m_clusterIndex.empty())
    {
      m_tree->Fill();
    }
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_evt
              << " poly_tracks=" << npoly_tracks
              << " residuals=" << nfilled << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_PolyClusterResiduals::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    m_outfile->cd();
    if (m_tree) m_tree->Write();
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
    m_tree = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
