#ifndef __TRACKJETINPUT_H__
#define __TRACKJETINPUT_H__

#include "JetInput.h"
#include "Jet.h"

#include <vector>

class PHCompositeNode;

class TrackJetInput : public JetInput {
  
public:

  TrackJetInput(Jet::SRC input);
  virtual ~TrackJetInput() {}

  void identify(std::ostream& os = std::cout);
  
  Jet::SRC get_src() {return _input;}
  
  std::vector<Jet*> get_input(PHCompositeNode *topNode);
  
private:
  int _verbosity;
  Jet::SRC _input;
};

#endif
