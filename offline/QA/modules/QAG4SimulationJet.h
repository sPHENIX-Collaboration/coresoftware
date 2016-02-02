#ifndef __CALOEVALUATOR_H__
#define __CALOEVALUATOR_H__

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <string>
#include <set>
#include <map>
#include <stdint.h>
#ifndef __CINT__
#include <memory>
#endif

#include <TString.h>

class PHCompositeNode;
class Fun4AllHistoManager;
class TH1F;
class JetEvalStack;

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
    kProcessComparison = 1 << 3,

    //! default. Do everything
    kDefaultFlag = kProcessTruthSpectrum | kProcessRecoSpectrum
        | kProcessComparison
  };

  QAG4SimulationJet(const std::string & truth_jet, enu_flags flags =
      kDefaultFlag);
  virtual
  ~QAG4SimulationJet();

  //! add reco jet to the process list
  //! @return number of reco jet on list
  int
  add_reco_jet(const std::string & reco_jet)
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

  int
  Init(PHCompositeNode *topNode);
  int
  InitRun(PHCompositeNode *topNode);
  int
  process_event(PHCompositeNode *topNode);
  int
  End(PHCompositeNode *topNode);

private:

  int
  Init_Spectrum(PHCompositeNode *topNode, const std::string & jet_name);
  int
  process_Spectrum(PHCompositeNode *topNode, const std::string & jet_name);

  int
  Init_Comparison(PHCompositeNode *topNode, const std::string & reco_jet_name);
  int
  process_Comparison(PHCompositeNode *topNode,
      const std::string & reco_jet_name);

  //! common prefix for QA histograms
  std::string
  get_histo_prefix(const std::string & src_jet_name = "",
      const std::string & reco_jet_name = "");

#ifndef __CINT__
  //! cache the jet evaluation modules
  typedef std::map<std::string, std::shared_ptr<JetEvalStack>> jetevalstacks_map;
  jetevalstacks_map _jetevalstacks;
#endif

  //! truth jet name
  std::string _truth_jet;

  //! list of reco jet
  std::set<std::string> _reco_jets;

  uint32_t _flags;

  //! simple counter
  unsigned long _ievent;

};

#endif // __CALOEVALUATOR_H__
