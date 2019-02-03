#ifndef __FASTJETALGOSUB_H__
#define __FASTJETALGOSUB_H__

#include <g4jets/Jet.h>
#include <g4jets/JetAlgo.h>

class FastJetAlgoSub : public JetAlgo {
  
public:

  FastJetAlgoSub(Jet::ALGO algo, float par, float verbosity = 0);
  virtual ~FastJetAlgoSub() {}

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
