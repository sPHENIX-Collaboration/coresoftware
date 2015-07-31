#ifndef __FASTJETALGO_H__
#define __FASTJETALGO_H__

#include <phool/PHCompositeNode.h>

#include "Jet.h"
#include "JetAlgo.h"

class FastJetAlgo : public JetAlgo {
  
public:

  FastJetAlgo();
  virtual ~FastJetAlgo() {}

  std::vector<Jet> get_jets(const std::vector<Jet>& particles);
  
protected:

private:
  int verbosity;
  
};

#endif
