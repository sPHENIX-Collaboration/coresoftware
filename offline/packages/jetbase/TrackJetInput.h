#ifndef G4JET_TRACKJETINPUT_H
#define G4JET_TRACKJETINPUT_H

#include "JetInput.h"

#include "Jet.h"

#include <iostream>  // for cout, ostream
#include <string>    // for string
#include <vector>

class PHCompositeNode;

class TrackJetInput : public JetInput
{
 public:
  TrackJetInput(Jet::SRC input, const std::string& name = "SvtxTrackMap");
  ~TrackJetInput() override {}

  void identify(std::ostream& os = std::cout) override;

  Jet::SRC get_src() override { return _input; }

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;

 private:
  std::string m_NodeName;
  Jet::SRC _input;
};

#endif
