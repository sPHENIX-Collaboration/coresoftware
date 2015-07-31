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

  virtual std::vector<Jet*> get_input(PHCompositeNode *topNode);

protected:

  
private:
  int verbosity;
};

#endif
