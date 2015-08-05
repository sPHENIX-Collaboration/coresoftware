#ifndef __TRUTHJETINPUT_H__
#define __TRUTHJETINPUT_H__

#include <phool/PHCompositeNode.h>

#include "JetInput.h"
#include "Jet.h"

#include <vector>

class TruthJetInput : public JetInput {
  
public:

  TruthJetInput();
  virtual ~TruthJetInput() {}

  std::vector<Jet*> get_input(PHCompositeNode *topNode);
  
  void set_eta_range(float eta_min, float eta_max) {
    _eta_min = eta_min;
    _eta_max = eta_max;
  }
    
private:
  int _verbosity;
  float _eta_min;
  float _eta_max;
};

#endif
