#ifndef G4EVAL_JETRECOEVAL_H
#define G4EVAL_JETRECOEVAL_H

#include "JetTruthEval.h"

#include <g4jets/Jet.h>

#include <map>
#include <set>
#include <string>
#include <utility>

class CaloEvalStack;

class JetMap;

class PHCompositeNode;

class PHG4Hit;
class PHG4Particle;
class PHG4Shower;

class RawClusterContainer;
class RawTowerContainer;

class SvtxEvalStack;
class SvtxTrackMap;

class JetRecoEval
{
 public:
  /// example recojetname:  AntiKt_Tower_r03
  /// example truthjetname: AntiKt_Truth_r03
  JetRecoEval(PHCompositeNode* topNode,
              const std::string& recojetname,
              const std::string& truthjetname);
  virtual ~JetRecoEval();

  /// reinitialize the eval for a new event
  void next_event(PHCompositeNode* topNode);

  /// activate or deactivate the memory caching inside the evaluation module
  void do_caching(bool do_cache)
  {
    _do_cache = do_cache;
    _jettrutheval.do_caching(do_cache);
  }

  /// strict mode will assert when an error is detected
  /// non-strict mode will notice and report at the End()
  void set_strict(bool strict)
  {
    _strict = strict;
    _jettrutheval.set_strict(strict);
  }

  /// get a count of the errors discovered thus far
  unsigned int get_errors() { return _errors + _jettrutheval.get_errors(); }

  /// adjust the messaging from the evalutaion module
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _jettrutheval.set_verbosity(verbosity);
  }

  /// get a copy of the lower level eval and its memory cache
  JetTruthEval* get_truth_eval() { return &_jettrutheval; }

  /// get a copy of the lower level eval and its memory cache
  SvtxEvalStack* get_svtx_eval_stack() { return _jettrutheval.get_svtx_eval_stack(); }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_cemc_eval_stack() { return _jettrutheval.get_cemc_eval_stack(); }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_hcalin_eval_stack() { return _jettrutheval.get_hcalin_eval_stack(); }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_hcalout_eval_stack() { return _jettrutheval.get_hcalout_eval_stack(); }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_femc_eval_stack() { return _jettrutheval.get_femc_eval_stack(); }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_fhcal_eval_stack() { return _jettrutheval.get_fhcal_eval_stack(); }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_eemc_eval_stack() { return _jettrutheval.get_eemc_eval_stack(); }
  // ---reduced sim node or better----------------------------------------------

  /// what truth showers contributed to this reconstructed jet?
  std::set<PHG4Shower*> all_truth_showers(Jet* recojet);

  /// what truth particles contributed to this reconstructed jet?
  std::set<PHG4Particle*> all_truth_particles(Jet* recojet);

  /// what truth jets contributed to this reconstructed jet?
  std::set<Jet*> all_truth_jets(Jet* recojet);

  /// which truth jet contributed the most energy to this reconstructed jet?
  Jet* max_truth_jet_by_energy(Jet* recojet);

  /// what reconstructed jets had contributions from this truth jet?
  std::set<Jet*> all_jets_from(Jet* truthjet);

  /// which reconstructed jet had the largest energy contribution from this truth jet?
  Jet* best_jet_from(Jet* truthjet);

  /// which reconstructed jet had the largest energy contribution from this truth jet in a unique match?
  /// @return pointer to reco jet. And NULL if no unique match
  Jet* unique_reco_jet_from_truth(Jet* truthjet);

  /// which truth jet had the largest energy contribution from this reco jet in a unique match?
  /// @return pointer to truth jet. And NULL if no unique match
  Jet* unique_truth_jet_from_reco(Jet* recojet);

  /// what was the energy contribution to this reconstructed jet from this truth jet?
  float get_energy_contribution(Jet* recojet, Jet* truthjet);

  /// what was the energy contribution to this reconstructed jet from a particular source
  float get_energy_contribution(Jet* recojet, Jet::SRC src);

  void set_track_nodename(const std::string& name);

  // ---full sim node required--------------------------------------------------

  /// which truth hits contributed to this reconstructed jet?
  std::set<PHG4Hit*> all_truth_hits(Jet* recojet);

 private:
  void get_node_pointers(PHCompositeNode* topNode);

  JetTruthEval _jettrutheval;
  std::string _recojetname;
  std::string _truthjetname;

  JetMap* _recojets;
  JetMap* _truthjets;

  SvtxTrackMap* _trackmap;
  RawTowerContainer* _cemctowers;
  RawClusterContainer* _cemcclusters;
  RawTowerContainer* _hcalintowers;
  RawClusterContainer* _hcalinclusters;
  RawTowerContainer* _hcalouttowers;
  RawClusterContainer* _hcaloutclusters;
  RawTowerContainer* _femctowers;
  RawClusterContainer* _femcclusters;
  RawTowerContainer* _fhcaltowers;
  RawClusterContainer* _fhcalclusters;
  RawTowerContainer* _eemctowers;
  RawClusterContainer* _eemcclusters;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::map<Jet*, std::set<PHG4Shower*> > _cache_all_truth_showers;
  std::map<Jet*, std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<Jet*, std::set<Jet*> > _cache_all_truth_jets;
  std::map<Jet*, Jet*> _cache_max_truth_jet_by_energy;
  std::map<Jet*, std::set<Jet*> > _cache_all_jets_from;
  std::map<Jet*, Jet*> _cache_best_jet_from;
  std::map<std::pair<Jet*, Jet*>, float> _cache_get_energy_contribution;
  std::map<std::pair<Jet*, Jet::SRC>, float> _cache_get_energy_contribution_src;  /// used in get_energy_contribution (Jet* recojet, Jet::SRC src);
  std::map<Jet*, std::set<PHG4Hit*> > _cache_all_truth_hits;
  std::string m_TrackNodeName = "SvtxTrackMap";
};

#endif  // G4EVAL_JETRECOEVAL_H
