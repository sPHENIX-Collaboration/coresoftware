#ifndef __FASTJETALGO_H__
#define __FASTJETALGO_H__

#include <phool/PHCompositeNode.h>

#include "Jet.h"
#include "JetAlgo.h"

class FastJetAlgo : public JetAlgo {
  
public:

  FastJetAlgo(Jet::ALGO algo, float par);
  virtual ~FastJetAlgo() {}

  Jet::ALGO get_algo() {return _algo;}
  float     get_par() {return _par;}
  
  std::vector<Jet*> get_jets(std::vector<Jet*> particles);
  
private:
  int _verbosity;
  Jet::ALGO _algo;
  float _par;
  
};

#endif
