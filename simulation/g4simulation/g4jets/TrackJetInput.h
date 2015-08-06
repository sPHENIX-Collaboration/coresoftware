#ifndef __TRACKJETINPUT_H__
#define __TRACKJETINPUT_H__

#include <phool/PHCompositeNode.h>

#include "JetInput.h"
#include "Jet.h"

#include <vector>

class TrackJetInput : public JetInput {
  
public:

  TrackJetInput(Jet::SRC input);
  virtual ~TrackJetInput() {}

  Jet::SRC get_src() {return _input;}
  
  std::vector<Jet*> get_input(PHCompositeNode *topNode);
  
private:
  int _verbosity;
  Jet::SRC _input;
};

#endif
