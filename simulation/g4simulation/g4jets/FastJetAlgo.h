#ifndef __FASTJETALGO_H__
#define __FASTJETALGO_H__

#include "Jet.h"
#include "JetAlgo.h"

// forward declarations
class PHCompositeNode;

class FastJetAlgo : public JetAlgo {
  
public:

  FastJetAlgo(Jet::ALGO algo, float par, float verbosity = 0);
  virtual ~FastJetAlgo() {}

  void      identify(std::ostream& os = std::cout);
  Jet::ALGO get_algo() {return _algo;}
  float     get_par() {return _par;}
  
  std::vector<Jet*> get_jets(std::vector<Jet*> particles);
  
private:
  int _verbosity;
  Jet::ALGO _algo;
  float _par;
  
};

#endif
