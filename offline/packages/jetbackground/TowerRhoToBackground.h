#ifndef JETBACKGROUND_TOWERRHOTOBACKGROUND_H
#define JETBACKGROUND_TOWERRHOTOBACKGROUND_H

#include <fun4all/SubsysReco.h>

#include <globalvertex/GlobalVertex.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;

/// \class TowerRhoToBackground
///
/// \brief Fills a TowerBackground from three TowerRho nodes (EMCal, IHCal, OHCal)
///
/// For each layer and eta strip the UE energy per tower is x = rho * cosh(eta), times
/// the strip delta eta * delta phi for the AREA method, where eta is the strip eta
/// corrected for the event vertex at the layer radius (as in TowerJetInput). The
/// EMCal layer is the retowered EMCal, on the IHCal grid at the EMCal radius.
/// For MULT, rho must come from the towers that are subtracted. v2 and Psi2 are 0.
/// The result is subtracted by SubtractTowers (set_inputTowerBackgroundNode).
///
/// Optionally the UE is calibrated per (layer, eta strip; vertex z, Q'), with
/// Q' = (q_S + q_N) / Q_knee(run): affine UE = alpha * x + beta (default) or
/// multiplicative UE = w * x. The calibration is a CDBHistos payload:
///   zvtx_bins, mbdq_bins  TH1, bin edges in z [cm] and Q', content = the
///                         event-weighted bin centre (the interpolation nodes)
///   mbdq_knee             TH1 over the run number, content = Q_knee (optional with a
///                         per-run knee payload, see set_mbd_knee_cdb_tag)
///   w_<L>, alpha_<L>, beta_<L>  TH3 (ieta, z, Q') for L = CEMC, HCALIN, HCALOUT
/// The calibration is interpolated bilinearly in (z, Q') between the bin centres,
/// clamped beyond the outermost centres (or looked up in the bin, see
/// set_interpolate). Events without a vertex or outside the calibrated vertex range
/// use UE = x, as do strips with an invalid calibration; all fallbacks are counted
/// (summary at End, getters below) and printed per event for Verbosity() > 1.
class TowerRhoToBackground : public SubsysReco
{
 public:
  enum class CalibModel
  {
    MULTIPLICATIVE,  // UE = w * x
    AFFINE           // UE = alpha * x + beta
  };

  explicit TowerRhoToBackground(const std::string &name = "TowerRhoToBackground")
    : SubsysReco(name)
  {
  }

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_emcal_rho_node(const std::string &name) { m_rho_nodes[0] = name; }
  void set_ihcal_rho_node(const std::string &name) { m_rho_nodes[1] = name; }
  void set_ohcal_rho_node(const std::string &name) { m_rho_nodes[2] = name; }
  void set_background_node(const std::string &name) { m_background_node = name; }

  void set_emcal_geom_node(const std::string &name) { m_geom_nodes[0] = name; }
  void set_ihcal_geom_node(const std::string &name) { m_geom_nodes[1] = name; }
  void set_ohcal_geom_node(const std::string &name) { m_geom_nodes[2] = name; }

  // UE calibration, off by default (UE = x): a CDB tag, or a direct path to the
  // calibration file (used until it is in the CDB); the path wins if both are set
  void set_eta_calib_cdb_tag(const std::string &tag) { m_eta_calib_tag = tag; }
  void set_eta_calib_path(const std::string &path) { m_eta_calib_path = path; }
  void set_calib_model(const CalibModel model) { m_model = model; }
  // bilinear interpolation between bin centres (default) or the value of the bin
  void set_interpolate(const bool val) { m_interpolate = val; }
  // MBD charge knee Q_knee of the run (MbdChargeKnee): set here, else from the per-run knee
  // payload (a CDB tag, or a path to one payload or to a directory of mbdq_knee_<run>.root; the
  // path wins), else from the calibration payload's per-run table, else the dataset default
  // (median of that table, with a warning)
  void set_mbd_knee(const float qknee) { m_knee_user = qknee; }
  void set_mbd_knee_cdb_tag(const std::string &tag) { m_knee_tag = tag; }
  void set_mbd_knee_path(const std::string &path) { m_knee_path = path; }
  // MBD node whose charge sum (both arms) gives Q
  void set_mbd_node(const std::string &name) { m_mbd_node = name; }

  // vertex for the eta correction and the calibration; by default the first vertex in the map
  void set_vertex_type(const GlobalVertex::VTXTYPE type)
  {
    m_vertex_type = type;
    m_use_vertex_type = true;
  }

  // counters
  uint64_t get_n_events() const { return m_n_events; }
  // events processed with no calibration configured: UE = x everywhere
  uint64_t get_n_events_no_calib() const { return m_n_no_calib; }
  // events without a vertex or outside the calibrated vertex range: UE = x everywhere
  uint64_t get_n_events_out_of_range() const { return m_n_out_of_range; }
  // events with a non-finite MBD charge sum (left uncalibrated)
  uint64_t get_n_events_bad_mbd_charge() const { return m_n_bad_mbdq; }
  // events beyond the outermost bin centres in z or Q' (calibration clamped, not a fallback)
  uint64_t get_n_events_zvtx_clamped() const { return m_n_zvtx_clamped; }
  uint64_t get_n_events_mbdq_clamped() const { return m_n_mbdq_clamped; }
  // (layer, eta strip) calibrations that were invalid (w <= 0 or NaN, alpha/beta NaN), and the events with any
  uint64_t get_n_invalid_weights() const { return m_n_invalid_w; }
  uint64_t get_n_events_invalid_weight() const { return m_n_events_invalid_w; }

 private:
  std::array<std::string, 3> m_rho_nodes{"TowerRho_CEMC_MULT", "TowerRho_HCALIN_MULT", "TowerRho_HCALOUT_MULT"};
  std::array<std::string, 3> m_geom_nodes{"TOWERGEOM_CEMC", "TOWERGEOM_HCALIN", "TOWERGEOM_HCALOUT"};
  std::string m_background_node{"TowerInfoBackground_Rho"};

  GlobalVertex::VTXTYPE m_vertex_type{GlobalVertex::UNDEFINED};
  bool m_use_vertex_type{false};

  // UE calibration, copied from the payload in InitRun
  std::string m_eta_calib_tag;
  std::string m_eta_calib_path;
  std::string m_mbd_node{"MbdOut"};
  CalibModel m_model{CalibModel::AFFINE};
  bool m_interpolate{true};
  float m_knee_user{-1};
  std::string m_knee_tag;
  std::string m_knee_path;
  float m_knee{-1};  // Q_knee in use for this run
  bool m_calib_on{false};
  int m_calib_neta{0};
  std::vector<float> m_zvtx_edges;
  std::vector<float> m_mbdq_edges;
  std::vector<float> m_zvtx_centres;
  std::vector<float> m_mbdq_centres;
  // per layer, index (izbin * n_mbdq + imbd) * n_eta + ieta
  std::array<std::vector<float>, 3> m_w{};
  std::array<std::vector<float>, 3> m_alpha{};
  std::array<std::vector<float>, 3> m_beta{};

  uint64_t m_n_events{0};
  uint64_t m_n_no_calib{0};
  uint64_t m_n_out_of_range{0};
  uint64_t m_n_bad_mbdq{0};
  uint64_t m_n_zvtx_clamped{0};
  uint64_t m_n_mbdq_clamped{0};
  uint64_t m_n_invalid_w{0};
  uint64_t m_n_events_invalid_w{0};

  int LoadEtaCalib();

  // interpolation nodes for val: indices i0, i1 into the centres and the weight of i1;
  // returns true if val is beyond the outermost centres (clamped)
  struct Nodes
  {
    int i0{0};
    int i1{0};
    float f{0};
  };
  static bool find_nodes(float val, const std::vector<float> &centres, const std::vector<float> &edges, bool interpolate, Nodes &nodes);
};

#endif
