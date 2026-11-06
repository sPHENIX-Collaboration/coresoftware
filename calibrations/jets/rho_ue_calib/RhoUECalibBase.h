#ifndef RHOUECALIB_RHOUECALIBBASE_H
#define RHOUECALIB_RHOUECALIBBASE_H

#include <globalvertex/GlobalVertex.h>

#include <fun4all/SubsysReco.h>

#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

class PHCompositeNode;
class TFile;

/// \class RhoUECalibBase
///
/// \brief Common part of the rho UE calibration makers (MbdChargeHistMaker,
/// TowerShapeMapMaker, RhoUECalibMaker)
///
/// Event selection: Au+Au minimum bias (MinimumBiasInfo, on by default) and a vertex
/// inside the class range in z. Event class = (vertex z bin, Q' bin) with
/// Q' = (q_S + q_N) / Q_knee(run); Q' outside the edges goes to the first / last bin.
/// The class binning defaults to the run 2024 Au+Au calibration (40 z x 26 Q').
/// The knee is set with set_mbd_knee() or read from a table (set_knee_table(), lines
/// "run knee ..." as written by RhoUECalibPayload::FitMbdKnees); makers that need no
/// classes (MbdChargeHistMaker) ignore it.
///
/// Every output histogram is a sum, so the files of all jobs of a dataset can be
/// merged with hadd. The knee used is stored as knee x events and events per run.
class RhoUECalibBase : public SubsysReco
{
 public:
  RhoUECalibBase(const std::string &name, const std::string &outfile);
  ~RhoUECalibBase() override;

  RhoUECalibBase(const RhoUECalibBase &) = delete;
  RhoUECalibBase &operator=(const RhoUECalibBase &) = delete;
  RhoUECalibBase(RhoUECalibBase &&) = delete;
  RhoUECalibBase &operator=(RhoUECalibBase &&) = delete;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // event selection
  void set_require_minbias(const bool val) { m_require_minbias = val; }
  void set_minbias_node(const std::string &name) { m_minbias_node = name; }
  void set_mbd_node(const std::string &name) { m_mbd_node = name; }
  // vertex for the class; by default the first vertex in the map
  void set_vertex_type(const GlobalVertex::VTXTYPE type)
  {
    m_vertex_type = type;
    m_use_vertex_type = true;
  }

  // classes
  void set_zvtx_edges(const std::vector<double> &edges) { m_zvtx_edges = edges; }
  void set_mbdq_edges(const std::vector<double> &edges) { m_mbdq_edges = edges; }
  void set_mbd_knee(const double qknee) { m_knee_user = qknee; }
  void set_knee_table(const std::string &file) { m_knee_table = file; }

  // towers: <prefix>_CEMC_RETOWER, <prefix>_HCALIN, <prefix>_HCALOUT
  void set_tower_prefix(const std::string &prefix) { m_tower_prefix = prefix; }
  // tower geometry (eta of the strips, layer radii), as in TowerRhoToBackground
  void set_geom_nodes(const std::string &cemc, const std::string &hcalin, const std::string &hcalout)
  {
    m_geom_nodes = {cemc, hcalin, hcalout};
  }

  static std::vector<double> default_zvtx_edges();
  static std::vector<double> default_mbdq_edges();

  static constexpr int NLAYER = 3;
  static constexpr int NETA = 24;
  static constexpr int NPHI = 64;
  static const std::array<std::string, NLAYER> &layer_names();  // CEMC, HCALIN, HCALOUT

 protected:
  // called for every selected event; zvtx, Q, Q' and the class are set
  virtual int fill(PHCompositeNode *topNode) = 0;
  // book the maker's histograms in the output file (current directory)
  virtual void book() = 0;
  // whether the maker needs the event class (and so the knee)
  virtual bool needs_class() const { return true; }

  int n_class() const { return n_zvtx() * n_mbdq(); }
  int n_zvtx() const { return static_cast<int>(m_zvtx_edges.size()) - 1; }
  int n_mbdq() const { return static_cast<int>(m_mbdq_edges.size()) - 1; }
  const std::vector<double> &zvtx_edges() const { return m_zvtx_edges; }
  const std::vector<double> &mbdq_edges() const { return m_mbdq_edges; }

  std::string tower_node(int layer) const;
  // whether the maker needs the vertex-corrected strip eta (m_eta_corr)
  virtual bool needs_geometry() const { return true; }

  // current selected event
  int run() const { return m_run; }
  int evt_index() const { return m_evt; }  // index of the selected event in this job (halves for precision checks)
  double zvtx() const { return m_zvtx; }
  double mbdq() const { return m_mbdq; }
  double mbdt0() const { return m_mbdt0; }
  double qprime() const { return m_qprime; }
  int izvtx() const { return m_izvtx; }
  int imbdq() const { return m_imbdq; }
  int event_class() const { return m_class; }  // izvtx * n_mbdq + imbdq
  // eta of strip ieta seen from the event vertex at the layer radius (layer 0 = the
  // retowered EMCal: IHCal strips at the EMCal radius), as in TowerRhoToBackground
  double eta_corr(const int layer, const int ieta) const { return m_eta_corr.at(layer).at(ieta); }

 private:
  int m_run{0};
  int m_evt{0};
  double m_zvtx{0};
  double m_mbdq{0};
  double m_mbdt0{0};
  double m_qprime{0};
  int m_izvtx{-1};
  int m_imbdq{-1};
  int m_class{-1};
  std::array<std::array<double, NETA>, NLAYER> m_eta_corr{};

  std::string m_outfile;
  TFile *m_file{nullptr};

  bool m_require_minbias{true};
  std::string m_minbias_node{"MinimumBiasInfo"};
  std::string m_mbd_node{"MbdOut"};
  GlobalVertex::VTXTYPE m_vertex_type{GlobalVertex::UNDEFINED};
  bool m_use_vertex_type{false};

  std::vector<double> m_zvtx_edges{default_zvtx_edges()};
  std::vector<double> m_mbdq_edges{default_mbdq_edges()};
  double m_knee_user{-1};
  std::string m_knee_table;
  std::map<int, double> m_knees;
  double m_knee{-1};

  std::string m_tower_prefix{"TOWERINFO_CALIB"};
  std::array<std::string, NLAYER> m_geom_nodes{"TOWERGEOM_CEMC", "TOWERGEOM_HCALIN", "TOWERGEOM_HCALOUT"};
  bool m_geom_loaded{false};
  std::array<double, NLAYER> m_radius{};
  std::array<std::array<double, NETA>, NLAYER> m_strip_eta{};
  int load_geometry(PHCompositeNode *topNode);

  uint64_t m_n_seen{0};
  uint64_t m_n_minbias{0};
  uint64_t m_n_vertex{0};
  uint64_t m_n_selected{0};
};

#endif
