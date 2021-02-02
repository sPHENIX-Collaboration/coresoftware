#ifndef G4JET_TRACKJETINPUT_H
#define G4JET_TRACKJETINPUT_H

#include "JetInput.h"

#include "Jet.h"

#include <iostream>    // for cout, ostream
#include <string>      // for string
#include <vector>

class PHCompositeNode;

class TrackJetInput : public JetInput
{
 public:
  TrackJetInput(Jet::SRC input, const std::string &name = "SvtxTrackMap");
  virtual ~TrackJetInput() {}

  void identify(std::ostream& os = std::cout);

  Jet::SRC get_src() { return _input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode);

 private:
  std::string m_NodeName;
  Jet::SRC _input;
};

#endif
