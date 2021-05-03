#ifndef QA_QAG4SIMULATIONJET_H
#define QA_QAG4SIMULATIONJET_H

#include <fun4all/SubsysReco.h>

#include <TString.h>

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>  // std::pair, std::make_pair

class JetEvalStack;
class JetTruthEval;
class Jet;
class PHCompositeNode;

/// \class QAG4SimulationJet
class QAG4SimulationJet : public SubsysReco
{
 public:
  enum enu_flags
  {
    //! spectrum of truth jets
    kProcessTruthSpectrum = 1 << 1,

    //! spectrum of reconstructed jets
    kProcessRecoSpectrum = 1 << 2,

    //! comparison of reco jet VS truth
    kProcessTruthMatching = 1 << 3,

    //! default. Do everything
    kDefaultFlag = kProcessTruthSpectrum | kProcessRecoSpectrum | kProcessTruthMatching
  };

  QAG4SimulationJet(const std::string &truth_jet, enu_flags flags =
                                                      kDefaultFlag);
  virtual ~QAG4SimulationJet() {}

  //! add reco jet to the process list
  //! @return number of reco jet on list
  int add_reco_jet(const std::string &reco_jet)
  {
    _reco_jets.insert(reco_jet);
    return _reco_jets.size();
  }

  uint32_t
  get_flags() const
  {
    return _flags;
  }

  void
  set_flags(enu_flags flags)
  {
    _flags = (uint32_t) flags;
  }

  void
  set_flag(enu_flags flag)
  {
    _flags |= (uint32_t) flag;
  }

  bool
  flag(enu_flags flag)
  {
    return _flags & flag;
  }

  void
  reset_flag(enu_flags flag)
  {
    _flags &= ~(uint32_t) flag;
  }

  //! Energy ratio difference cut from 1 for matched jets
  double
  get_jet_match_dE_Ratio() const
  {
    return _jet_match_dE_Ratio;
  }

  //! Energy ratio difference cut from 1 for matched jets
  void
  set_jet_match_dE_Ratio(double jetMatchDERatio)
  {
    _jet_match_dE_Ratio = jetMatchDERatio;
  }

  //! Eta difference cut for matched jets
  double
  get_jet_match_dEta() const
  {
    return _jet_match_dEta;
  }

  //! Eta difference cut for matched jets
  void
  set_jet_match_dEta(double jetMatchDEta)
  {
    _jet_match_dEta = jetMatchDEta;
  }

  //! Phi difference cut for matched jets
  double
  get_jet_match_dPhi() const
  {
    return _jet_match_dPhi;
  }

  //! Phi difference cut for matched jets
  void
  set_jet_match_dPhi(double jetMatchDPhi)
  {
    _jet_match_dPhi = jetMatchDPhi;
  }

  //! set eta range
  void
  set_eta_range(double low, double high);

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

 private:
  int Init_Spectrum(PHCompositeNode *topNode, const std::string &jet_name);
  int process_Spectrum(PHCompositeNode *topNode, const std::string &jet_name, const bool is_reco_jet);

  int Init_TruthMatching(PHCompositeNode *topNode, const std::string &reco_jet_name);
  int process_TruthMatching(PHCompositeNode *topNode,
                            const std::string &reco_jet_name);

  //! common prefix for QA histograms
  std::string
  get_histo_prefix(const std::string &src_jet_name = "",
                   const std::string &reco_jet_name = "");

  //! cache the jet evaluation modules
  typedef std::map<std::string, std::shared_ptr<JetEvalStack>> jetevalstacks_map;
  jetevalstacks_map _jetevalstacks;
  std::shared_ptr<JetTruthEval> _jettrutheval;

  //! truth jet name
  std::string _truth_jet;

  //! list of reco jet
  std::set<std::string> _reco_jets;

  uint32_t _flags;

  //! eta range
  std::pair<double, double> eta_range;

  //! string description of eta range
  //! @return TString as ROOT likes
  TString
  get_eta_range_str(const char *eta_name = "#eta_{Jet}") const;

  //! acceptance cut on jet object
  bool
  jet_acceptance_cut(const Jet *jet) const;

  //! Eta difference cut for matched jets
  double _jet_match_dEta;

  //! Phi difference cut for matched jets
  double _jet_match_dPhi;

  //! Energy ratio difference cut from 1 for matched jets
  double _jet_match_dE_Ratio;
};

#endif  // QA_QAG4SIMULATIONJET_H
