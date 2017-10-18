#ifndef __TOWERJETINPUT_H__
#define __TOWERJETINPUT_H__

#include <phool/PHCompositeNode.h>

#include "JetInput.h"
#include "Jet.h"

#include <vector>

class TowerJetInput : public JetInput {
  
public:

  TowerJetInput(Jet::SRC input, double ecut, double towerscale);
  virtual ~TowerJetInput() {}

  void identify(std::ostream& os = std::cout);
  
  Jet::SRC get_src() {return _input;}
  
  std::vector<Jet*> get_input(PHCompositeNode *topNode);
  
private:
  int _verbosity;
  Jet::SRC _input;
  double _ecut;
  double _towerscale;
};

#endif
