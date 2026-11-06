#ifndef RHOUECALIB_RHOUECALIBEVENTVETO_H
#define RHOUECALIB_RHOUECALIBEVENTVETO_H

#include <globalvertex/GlobalVertex.h>

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <map>
#include <string>

class PHCompositeNode;

/// \class RhoUECalibEventVeto
///
/// \brief Event veto in front of the rho UE calibration makers, or of a DST skim (ABORTEVENT
/// for vetoed events): Au+Au minimum bias (optional), a vertex of the chosen type with
/// |z| < z_max, the bunch crossing, and no MBD charge loss (off-band pile-up). Every cut can
/// be switched off.
///
/// Bunch: MinimumBiasInfo::getBunchClass() == MAIN (or a satellite with
/// set_accept_satellites) when the MinimumBiasClassifier ran with bunch classes
/// (getT0Corrected() finite); otherwise |t0' = t0 - t0_offset(run)| < t0_window with the MbdOut
/// t0 and the run's main-bunch centre from the knee table (column "t0", written by
/// RhoUECalibPayload::FitMbdKnees) or set_t0_offset(). With set_tag_bunch_class the t0' and
/// its class (MAIN |t0'| < main window, M5 / P5 up to the satellite window, else NONE) are
/// written into MinimumBiasInfo (v2) when the classifier left them unset.
/// Off-band: full MBD charge (sum of the PMT charges > 0, MbdPmtContainer, reconstructed from
/// the waveforms by the current MbdReco, which keeps the charge of PMTs without a time hit) /
/// the production MBD charge > offband_cut. The production charge is the MbdOut of the
/// original DST_CALOFITTING, read into a second top node (set_production_topnode, e.g. a
/// Fun4AllDstInputManager("...", "DST", "TOPPROD") with BranchSelect of DST#MBD#MbdOut);
/// without it the refit MbdOut is used, which already contains most PMTs without a time hit
/// (charge-channel times), so the cut rejects almost nothing.
class RhoUECalibEventVeto : public SubsysReco
{
 public:
  explicit RhoUECalibEventVeto(const std::string &name = "RhoUECalibEventVeto");
  ~RhoUECalibEventVeto() override = default;
  RhoUECalibEventVeto(const RhoUECalibEventVeto &) = delete;
  RhoUECalibEventVeto &operator=(const RhoUECalibEventVeto &) = delete;
  RhoUECalibEventVeto(RhoUECalibEventVeto &&) = delete;
  RhoUECalibEventVeto &operator=(RhoUECalibEventVeto &&) = delete;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // Au+Au minimum bias (MinimumBiasInfo::isAuAuMinimumBias), off by default
  void set_require_minbias(const bool val) { m_require_minbias = val; }

  // vertex: a vertex of this type (default: the first in the map) with |z| < z_max (< 0: off)
  void set_vertex_type(const GlobalVertex::VTXTYPE type)
  {
    m_vertex_type = type;
    m_use_vertex_type = true;
  }
  void set_zvtx_max(const double zmax) { m_zvtx_max = zmax; }

  // main bunch (off by default): t0 window around the run's main-bunch centre
  void set_main_bunch(const bool val) { m_main_bunch = val; }
  void set_t0_window(const double ns) { m_t0_window = ns; }
  void set_t0_offset(const double ns) { m_t0_offset_user = ns; }
  void set_t0_table(const std::string &file) { m_t0_table = file; }  // knee table with a "t0" column
  // classifier bunch classes: also accept the +-5 ns satellites (table mode: use the window)
  void set_accept_satellites(const bool val) { m_accept_satellites = val; }
  // write t0' and the bunch class into MinimumBiasInfo when the classifier left them unset
  void set_tag_bunch_class(const bool val) { m_tag_bunch_class = val; }
  void set_bunch_class_windows(const double main_ns, const double satellite_ns)
  {
    m_main_window = main_ns;
    m_satellite_window = satellite_ns;
  }

  // off-band: full refit / production charge > cut (<= 0: off)
  void set_offband_cut(const double cut) { m_offband_cut = cut; }
  // top node holding the production MbdOut (same node name as set_mbd_node)
  void set_production_topnode(const std::string &name) { m_prod_topnode = name; }

  void set_minbias_node(const std::string &name) { m_minbias_node = name; }
  void set_mbd_node(const std::string &name) { m_mbd_node = name; }
  void set_mbd_pmt_node(const std::string &name) { m_mbd_pmt_node = name; }

 private:
  GlobalVertex::VTXTYPE m_vertex_type{GlobalVertex::UNDEFINED};
  bool m_use_vertex_type{false};
  double m_zvtx_max{60.};

  bool m_require_minbias{false};
  bool m_main_bunch{false};
  double m_t0_window{2.5};
  bool m_accept_satellites{false};
  bool m_tag_bunch_class{false};
  double m_main_window{2.5};
  double m_satellite_window{7.5};
  double m_t0_offset_user{0.};
  std::string m_t0_table;
  std::map<int, double> m_t0_offsets;
  double m_t0_offset{0.};
  bool m_t0_offset_from_table{false};

  double m_offband_cut{1.3};
  std::string m_prod_topnode;

  std::string m_minbias_node{"MinimumBiasInfo"};
  std::string m_mbd_node{"MbdOut"};
  std::string m_mbd_pmt_node{"MbdPmtContainer"};

  uint64_t m_n_seen{0};
  uint64_t m_n_minbias{0};
  uint64_t m_n_vertex{0};
  uint64_t m_n_bunch{0};
  uint64_t m_n_bunch_classifier{0};
  uint64_t m_n_offband_veto{0};
  uint64_t m_n_no_prod_charge{0};
  uint64_t m_n_selected{0};
};

#endif
