#ifndef G4EVAL_JETTRUTHEVAL_H
#define G4EVAL_JETTRUTHEVAL_H

#include "CaloEvalStack.h"
#include "SvtxEvalStack.h"

#include <map>
#include <set>
#include <string>

class Jet;
class JetMap;

class PHCompositeNode;

class PHG4Hit;
class PHG4Particle;
class PHG4Shower;
class PHG4TruthInfoContainer;

class JetTruthEval
{
 public:
  /// example truthjetname: AntiKt_Truth_r03
  JetTruthEval(PHCompositeNode* topNode,
               const std::string& truthjetname);
  virtual ~JetTruthEval();

  /// reinitialize the eval for a new event
  void next_event(PHCompositeNode* topNode);

  /// activate or deactivate the memory caching inside the evaluation module
  void do_caching(bool do_cache)
  {
    _do_cache = do_cache;
    _svtxevalstack.do_caching(do_cache);
    _cemcevalstack.do_caching(do_cache);
    _hcalinevalstack.do_caching(do_cache);
    _hcaloutevalstack.do_caching(do_cache);
    _femcevalstack.do_caching(do_cache);
    _fhcalevalstack.do_caching(do_cache);
    _eemcevalstack.do_caching(do_cache);
  }

  /// strict mode will assert when an error is detected
  /// non-strict mode will notice and report at the End()
  void set_strict(bool strict)
  {
    _strict = strict;
    _svtxevalstack.set_strict(strict);
    _cemcevalstack.set_strict(strict);
    _hcalinevalstack.set_strict(strict);
    _hcaloutevalstack.set_strict(strict);
    _femcevalstack.set_strict(strict);
    _fhcalevalstack.set_strict(strict);
    _eemcevalstack.set_strict(strict);
  }

  /// get a count of the errors discovered thus far
  unsigned int get_errors()
  {
    return _errors + _svtxevalstack.get_errors() + _cemcevalstack.get_errors() + _hcalinevalstack.get_errors() + _hcaloutevalstack.get_errors() + _femcevalstack.get_errors() + _fhcalevalstack.get_errors() + _eemcevalstack.get_errors();
  }

  /// adjust the messaging from the evalutaion module
  void set_verbosity(int verbosity)
  {
    _verbosity = verbosity;
    _svtxevalstack.set_verbosity(verbosity);
    _cemcevalstack.set_verbosity(verbosity);
    _hcalinevalstack.set_verbosity(verbosity);
    _hcaloutevalstack.set_verbosity(verbosity);
    _femcevalstack.set_verbosity(verbosity);
    _fhcalevalstack.set_verbosity(verbosity);
    _eemcevalstack.set_verbosity(verbosity);
  }

  //set track node name
  void set_track_nodename(const std::string& name)
  {
    _svtxevalstack.set_track_nodename(name);
  }

  /// get a copy of the lower level eval and its memory cache
  SvtxEvalStack* get_svtx_eval_stack() { return &_svtxevalstack; }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_cemc_eval_stack() { return &_cemcevalstack; }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_hcalin_eval_stack() { return &_hcalinevalstack; }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_hcalout_eval_stack() { return &_hcaloutevalstack; }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_femc_eval_stack() { return &_femcevalstack; }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_fhcal_eval_stack() { return &_fhcalevalstack; }

  /// get a copy of the lower level eval and its memory cache
  CaloEvalStack* get_eemc_eval_stack() { return &_eemcevalstack; }

  // ---reduced sim node or better----------------------------------------------

  /// which truth jet in the specified node contains this truth particle?
  Jet* get_truth_jet(PHG4Particle* truthparticle);

  /// which truth particle contributed to this truth jet?
  std::set<PHG4Particle*> all_truth_particles(Jet* truthjet);

  /// which showers were left by particles contributing to this truth jet?
  std::set<PHG4Shower*> all_truth_showers(Jet* truthjet);

  // ---full sim node required--------------------------------------------------

  /// which truth hits were left by particles contributing to this truth jet?
  std::set<PHG4Hit*> all_truth_hits(Jet* truthjet);

 private:
  void get_node_pointers(PHCompositeNode* topNode);

  std::string _truthjetname;
  SvtxEvalStack _svtxevalstack;
  CaloEvalStack _cemcevalstack;
  CaloEvalStack _hcalinevalstack;
  CaloEvalStack _hcaloutevalstack;
  CaloEvalStack _femcevalstack;
  CaloEvalStack _fhcalevalstack;
  CaloEvalStack _eemcevalstack;
  
  PHG4TruthInfoContainer* _truthinfo;
  JetMap* _truthjets;

  bool _strict;
  int _verbosity;
  unsigned int _errors;

  bool _do_cache;
  std::map<Jet*, std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<Jet*, std::set<PHG4Shower*> > _cache_all_truth_showers;
  std::map<Jet*, std::set<PHG4Hit*> > _cache_all_truth_hits;
  std::map<PHG4Particle*, Jet*> _cache_get_truth_jet;
};

#endif  // JETTRUTHEVAL_H__
