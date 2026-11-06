/*!
 * \file DetermineTowerRho.h
 * \brief UE background rho calculator.
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 2.0.1 $
 * \date $Date: 02/01/2024. Revised 09/19/2024$
 */

#ifndef JETBACKGROUND_DETERMINETOWERHO_H
#define JETBACKGROUND_DETERMINETOWERHO_H

#include "TowerRho.h"

#include <globalvertex/GlobalVertex.h>

#include <fun4all/SubsysReco.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

class PHCompositeNode;
class JetContainer;
class JetInput;
class Jet;

namespace fastjet
{
  class PseudoJet;
  class Selector;
}  // namespace fastjet

class DetermineTowerRho : public SubsysReco
{
 public:
  DetermineTowerRho(const std::string &name = "DetermineTowerRho");
  ~DetermineTowerRho() override;

  // standard Fun4All methods
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  // add rho method (Area or Multiplicity)
  //
  // The background jets are clustered with fastjet and selected with the jet
  // acceptance (set_jet_abs_eta). The surviving jets are then converted back
  // to sPHENIX Jets: the four-vector of each is rebuilt from the ORIGINAL
  // constituents,
  // Every selected jet -- the omitted hardest ones (the seeds) included -- is
  // written to the node tree as a JetContainer named jet_node. The default rho
  // node is TowerRho_<METHOD> and the default jet node is <rho node>_SeedJets,
  // e.g. TowerRho_MULT -> TowerRho_MULT_SeedJets, or with output_node =
  // "TowerRho_CEMC_MULT" -> TowerRho_CEMC_MULT_SeedJets.
  //   Jet::PROPERTY::prop_SeedItr : 1 if it was omitted as one of the hardest
  //                                 (a seed), 0 if it entered rho
  //   Jet::PROPERTY::prop_signedSumeT : the signed scalar sum of its constituent
  //                                 pT, the quantity rho is estimated from

  void add_method(
    TowerRho::Method rho_method,
    const std::string &output_node = "",
    const std::string &jet_node = ""
  );

  // inputs for estimating background
  void add_input(JetInput *input) { m_inputs.push_back(input); }
  void add_tower_input(JetInput *input) { add_input(input); }  // for backwards compatibility

  // -- background jet clustering --
  // algorithm (default KT) and radius parameter R (default 0.4)
  void set_algo(const Jet::ALGO algo) { m_bkgd_jet_algo = algo; }
  Jet::ALGO get_algo() const { return m_bkgd_jet_algo; }

  void set_par(const float val) { m_par = val; }
  float get_par() const { return m_par; }

  // -- ghosts (AREA method only) --
  // ghost area (default 0.01)
  void set_ghost_area(const float val) { m_ghost_area = val; }
  float get_ghost_area() const { return m_ghost_area; }

  // -- background jet acceptance --
  void set_jet_abs_eta(const float abseta) { m_abs_jet_eta_range = abseta; }
  float get_jet_abs_eta() const { return m_abs_jet_eta_range; }

  // number of hardest jets, ranked by prop_signedSumeT, left out of the rho
  // estimate (default 2)
  void set_omit_nhardest(const unsigned int val) { m_omit_nhardest = val; }
  unsigned int get_omit_nhardest() const { return m_omit_nhardest; }

  // -- tower shape map (s-map), off by default --
  // s(class, layer, ieta, iphi) = <E_T,tower> / <E_T> over the towers of the layer, per class of
  // vertex z and Q' = (q_S + q_N) / Q_knee(run). With a map set, every tower is divided by its s
  // for the background-jet clustering and for the choice of the omitted hardest jets (the
  // seeds); rho and sigma are still taken from the unscaled towers ("seed only"), and the jets
  // written to the jet node keep their unscaled four-vectors and prop_signedSumeT.
  // Works with TOWERINFO inputs of the retowered EMCal, the IHCal and the OHCal (24 x 64).
  // Events without a vertex, or outside the map's vertex range, use s = 1 (counted).
  // CDBHistos payload: zvtx_bins and mbdq_bins (TH1, bin edges; Q' edges), mbdq_knee (TH1 over the
  // run number, content Q_knee; optional with a knee payload), s_CEMC, s_HCALIN, s_HCALOUT (TH3:
  // ieta 24, iphi 64, class = izbin * n_mbdQ_bins + imbd); s <= 0 or NaN means s = 1. A CDB tag,
  // or a direct path to the payload (used until it is in the CDB); the path wins if both are set.
  void set_tower_shape_cdb_tag(const std::string &tag) { m_smap_tag = tag; }
  void set_tower_shape_path(const std::string &path) { m_smap_path = path; }
  // MBD charge knee Q_knee of the run (MbdChargeKnee): set here, else from the per-run knee
  // payload (a CDB tag, or a path to one payload or to a directory of mbdq_knee_<run>.root; the
  // path wins), else from the s-map payload's per-run table, else the dataset default (median of
  // that table, with a warning)
  void set_mbd_knee(const float qknee) { m_smap_knee_user = qknee; }
  void set_mbd_knee_cdb_tag(const std::string &tag) { m_knee_tag = tag; }
  void set_mbd_knee_path(const std::string &path) { m_knee_path = path; }
  // MBD node whose charge sum (both arms) gives Q
  void set_mbd_node(const std::string &name) { m_mbd_node = name; }
  // vertex for the map class; by default the first vertex in the map
  void set_vertex_type(const GlobalVertex::VTXTYPE type)
  {
    m_vertex_type = type;
    m_use_vertex_type = true;
  }

  // s-map counters: events processed with the map, and those that used s = 1
  uint64_t get_n_events_tower_shape() const { return m_n_smap_events; }
  uint64_t get_n_events_tower_shape_no_vertex() const { return m_n_smap_novtx; }
  uint64_t get_n_events_tower_shape_out_of_range() const { return m_n_smap_out_of_range; }

 private:

  // variables
  std::vector<JetInput *> m_inputs;

  std::vector<std::string> m_output_nodes;
  std::vector<std::string> m_jet_output_nodes;

  std::vector<TowerRho::Method> m_rho_methods;

  Jet::ALGO m_bkgd_jet_algo{Jet::ALGO::KT};
  float m_par{0.4};

  float m_ghost_area{0.01};

  float m_abs_jet_eta_range{-1.0};  // negative = event tower eta range shrunk by R
  unsigned int m_omit_nhardest{2};

  // tower shape map
  std::string m_smap_tag;
  std::string m_smap_path;
  std::string m_mbd_node{"MbdOut"};
  GlobalVertex::VTXTYPE m_vertex_type{GlobalVertex::UNDEFINED};
  bool m_use_vertex_type{false};
  float m_smap_knee_user{-1};
  std::string m_knee_tag;
  std::string m_knee_path;
  float m_smap_knee{-1};  // Q_knee in use for this run
  bool m_smap_on{false};
  int m_smap_neta{0};
  int m_smap_nphi{0};
  std::vector<float> m_smap_zvtx_edges;
  std::vector<float> m_smap_mbdq_edges;
  // per layer (0 retowered EMCal, 1 IHCal, 2 OHCal), indexed like the payload channel
  std::array<std::vector<float>, 3> m_smap{};
  uint64_t m_n_smap_events{0};
  uint64_t m_n_smap_novtx{0};
  uint64_t m_n_smap_out_of_range{0};

  // internal methods
  int CreateNodes(PHCompositeNode *topNode);

  // load the s-map for this run (InitRun); ABORTRUN if it is set but unusable
  int LoadTowerShape();
  // map class of this event, -1 (s = 1) without a vertex or out of range
  int TowerShapeClass(PHCompositeNode *topNode);
  // s of every input particle for map class iclass (1 if iclass < 0)
  std::vector<float> TowerShapeScales(const std::vector<Jet *> &particles, int iclass) const;
  // layer of a TOWERINFO input in the map, -1 if it has none
  static int TowerShapeLayer(int src);

  // Convert the selected fastjet jets back to sPHENIX Jets; rank_pt receives, per jet, the
  // signed sum of the constituent pT divided by their s (= prop_signedSumeT without a map)
  void ConvertJets(
    JetContainer *jets,
    const std::vector<fastjet::PseudoJet> &fastjets,
    const std::vector<Jet *> &particles,
    const std::vector<float> &scales,
    bool with_area,
    std::vector<float> &rank_pt
  ) const;


  // prop_signedSumeT of jet ijet: the per-jet quantity the median is taken over,
  static float SignedSumeT(JetContainer *jets, unsigned int ijet);

  // Drop the m_omit_nhardest hardest jets, ranked by rank_pt
  std::vector<unsigned int> SelectJets(const std::vector<float> &rank_pt) const;

  // Estimate rho and sigma from the sPHENIX jets of `jets` listed in `keep`
  void CalcRho(
    JetContainer *jets,
    const std::vector<unsigned int> &keep,
    TowerRho::Method rho_method,
    float &rho,
    float &sigma
  ) const;

  static float CalcPercentile(
    const std::vector<float> &sorted_vec,
    float percentile,
    float nempty
  );

  static void CalcMedianStd(
    const std::vector<float> &vec,
    float n_empty_jets,
    float &median,
    float &std_dev
  );

  fastjet::Selector get_jet_selector(double eta_min, double eta_max) const;
};

#endif  // JETBACKGROUND_DETERMINETOWERHO_H
