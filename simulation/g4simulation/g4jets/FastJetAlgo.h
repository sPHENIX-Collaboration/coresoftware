#ifndef __FASTJETALGO_H__
#define __FASTJETALGO_H__

#include <phool/PHCompositeNode.h>

#include "Jet.h"
#include "JetAlgo.h"

class FastJetAlgo : public JetAlgo {
  
public:

  enum {ANTIKT=0,KT=1,CA=2};
  
  FastJetAlgo(int algo, float par);
  virtual ~FastJetAlgo() {}

  std::vector<Jet*> get_jets(std::vector<Jet*> particles);
  
private:
  int _verbosity;
  int _algo;
  float _par;
  
};

#endif
