
#ifndef __JETTRUTHEVAL_H__
#define __JETTRUTHEVAL_H__

#include "SvtxEvalStack.h"
#include "CaloEvalStack.h"

#include <phool/PHCompositeNode.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4Shower.h>
#include <g4main/PHG4Hit.h>
#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

#include <map>
#include <string>

class JetTruthEval {

public:

  JetTruthEval(PHCompositeNode* topNode,
	       std::string truthjetname);
  virtual ~JetTruthEval();

  void next_event(PHCompositeNode* topNode);
  void do_caching(bool do_cache) {
    _do_cache = do_cache;
    _svtxevalstack.do_caching(do_cache);
    _cemcevalstack.do_caching(do_cache);
    _hcalinevalstack.do_caching(do_cache);
    _hcaloutevalstack.do_caching(do_cache);
    _femcevalstack.do_caching(do_cache);
    _fhcalevalstack.do_caching(do_cache);
  }
  void set_strict(bool strict) {
    _strict = strict;
    _svtxevalstack.set_strict(strict);
    _cemcevalstack.set_strict(strict);
    _hcalinevalstack.set_strict(strict);
    _hcaloutevalstack.set_strict(strict); 
    _femcevalstack.set_strict(strict); 
    _fhcalevalstack.set_strict(strict); 
  }  
  void set_verbosity(int verbosity) {
    _verbosity = verbosity;
    _svtxevalstack.set_verbosity(verbosity);
    _cemcevalstack.set_verbosity(verbosity);
    _hcalinevalstack.set_verbosity(verbosity);
    _hcaloutevalstack.set_verbosity(verbosity); 
    _femcevalstack.set_verbosity(verbosity); 
    _fhcalevalstack.set_verbosity(verbosity); 
  }
  
  SvtxEvalStack* get_svtx_eval_stack() {return &_svtxevalstack;}
  CaloEvalStack* get_cemc_eval_stack() {return &_cemcevalstack;}
  CaloEvalStack* get_hcalin_eval_stack() {return &_hcalinevalstack;}
  CaloEvalStack* get_hcalout_eval_stack() {return &_hcaloutevalstack;}
  CaloEvalStack* get_femc_eval_stack() {return &_femcevalstack;}
  CaloEvalStack* get_fhcal_eval_stack() {return &_fhcalevalstack;}

  std::set<PHG4Particle*> all_truth_particles(Jet* truthjet);
  std::set<PHG4Shower*>   all_truth_showers(Jet* truthjet);  
  std::set<PHG4Hit*>      all_truth_hits(Jet* truthjet);
  
  Jet* get_truth_jet(PHG4Particle* truthparticle);
  
  unsigned int get_errors() {
    return _errors
      + _svtxevalstack.get_errors()
      + _cemcevalstack.get_errors()
      + _hcalinevalstack.get_errors()
      + _hcaloutevalstack.get_errors()
      + _femcevalstack.get_errors()
      + _fhcalevalstack.get_errors();
  }
  
private:

  void get_node_pointers(PHCompositeNode* topNode);

  std::string _truthjetname;
  SvtxEvalStack _svtxevalstack;
  CaloEvalStack _cemcevalstack;
  CaloEvalStack _hcalinevalstack;
  CaloEvalStack _hcaloutevalstack;
  CaloEvalStack _femcevalstack;
  CaloEvalStack _fhcalevalstack;

  PHG4TruthInfoContainer* _truthinfo;
  JetMap* _truthjets;

  bool _strict;
  int _verbosity;
  unsigned int _errors;
  
  bool                                    _do_cache;
  std::map<Jet*,std::set<PHG4Particle*> > _cache_all_truth_particles;
  std::map<Jet*,std::set<PHG4Shower*> >   _cache_all_truth_showers;
  std::map<Jet*,std::set<PHG4Hit*> >      _cache_all_truth_hits;
  std::map<PHG4Particle*,Jet*>            _cache_get_truth_jet;  
};

#endif // __JETTRUTHEVAL_H__
