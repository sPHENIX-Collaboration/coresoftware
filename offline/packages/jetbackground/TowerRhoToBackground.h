#ifndef JETBACKGROUND_TOWERRHOTOBACKGROUND_H
#define JETBACKGROUND_TOWERRHOTOBACKGROUND_H

#include <fun4all/SubsysReco.h>

#include <globalvertex/GlobalVertex.h>

#include <array>
#include <string>
#include <vector>

class CDBTTree;
class PHCompositeNode;

/// \class TowerRhoToBackground
///
/// \brief Fills a TowerBackground from three TowerRho nodes (EMCal, IHCal, OHCal)
///
/// For each layer and eta strip the UE energy per tower is rho * cosh(eta), times
/// the strip delta eta * delta phi for the AREA method, where eta is the strip eta
/// corrected for the event vertex at the layer radius (as in TowerJetInput). The
/// EMCal layer is the retowered EMCal, on the IHCal grid at the EMCal radius.
/// For MULT, rho must come from the towers that are subtracted. v2 and Psi2 are 0.
/// The result is subtracted by SubtractTowers (set_inputTowerBackgroundNode).
///
/// Optionally the UE is multiplied by an eta-shape calibration
/// w(layer, ieta; vertex z, MBD charge), read from a CDBTTree with single values
/// n_eta, n_zvtx_bins, n_mbdQ_bins, zvtx_edge_<i>, mbdQ_edge_<i> and per-channel
/// fields w_cemc / w_hcalin / w_hcalout, channel = izbin * (n_mbdQ_bins * n_eta) +
/// imbd * n_eta + ieta. Events outside the calibrated vertex z / MBD charge range,
/// and non-positive weights, use w = 1.
class TowerRhoToBackground : public SubsysReco
{
 public:
  TowerRhoToBackground(const std::string &name = "TowerRhoToBackground")
    : SubsysReco(name)
  {
  }
  ~TowerRhoToBackground() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void set_emcal_rho_node(const std::string &name) { m_rho_nodes[0] = name; }
  void set_ihcal_rho_node(const std::string &name) { m_rho_nodes[1] = name; }
  void set_ohcal_rho_node(const std::string &name) { m_rho_nodes[2] = name; }
  void set_background_node(const std::string &name) { m_background_node = name; }

  void set_emcal_geom_node(const std::string &name) { m_geom_nodes[0] = name; }
  void set_ihcal_geom_node(const std::string &name) { m_geom_nodes[1] = name; }
  void set_ohcal_geom_node(const std::string &name) { m_geom_nodes[2] = name; }

  // eta-shape calibration, off by default (w = 1): a CDB tag, or a direct path to
  // the calibration file (used until it is in the CDB); the path wins if both are set
  void set_eta_calib_cdb_tag(const std::string &tag) { m_eta_calib_tag = tag; }
  void set_eta_calib_path(const std::string &path) { m_eta_calib_path = path; }
  // MBD node whose charge sum (both arms) indexes the calibration
  void set_mbd_node(const std::string &name) { m_mbd_node = name; }

  // vertex for the eta correction; by default the first vertex in the map
  void set_vertex_type(const GlobalVertex::VTXTYPE type)
  {
    m_vertex_type = type;
    m_use_vertex_type = true;
  }

 private:
  std::array<std::string, 3> m_rho_nodes{"TowerRho_CEMC_MULT", "TowerRho_HCALIN_MULT", "TowerRho_HCALOUT_MULT"};
  std::array<std::string, 3> m_geom_nodes{"TOWERGEOM_CEMC", "TOWERGEOM_HCALIN", "TOWERGEOM_HCALOUT"};
  std::string m_background_node{"TowerInfoBackground_Rho"};

  GlobalVertex::VTXTYPE m_vertex_type{GlobalVertex::UNDEFINED};
  bool m_use_vertex_type{false};

  // eta-shape calibration
  std::string m_eta_calib_tag{};
  std::string m_eta_calib_path{};
  std::string m_mbd_node{"MbdOut"};
  CDBTTree *m_eta_calib{nullptr};
  int m_calib_neta{0};
  std::vector<float> m_calib_zvtx_edges{};
  std::vector<float> m_calib_mbdq_edges{};

  int LoadEtaCalib();
  // index of the bin of val in edges, -1 if outside
  static int find_bin(float val, const std::vector<float> &edges);
};

#endif
