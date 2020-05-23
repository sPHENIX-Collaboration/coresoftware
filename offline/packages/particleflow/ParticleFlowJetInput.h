#ifndef G4JET_PARTICLEFLOWJETINPUT_H
#define G4JET_PARTICLEFLOWJETINPUT_H

// finally system includes
#include <iostream>    // for cout, ostream
#include <vector>

#include <g4jets/Jet.h>
#include <g4jets/JetInput.h>

// forward declarations
class PHCompositeNode;

class ParticleFlowJetInput : public JetInput
{
 public:
  ParticleFlowJetInput();
  virtual ~ParticleFlowJetInput() {}

  std::vector<Jet*> get_input(PHCompositeNode* topNode);
  void identify(std::ostream& os = std::cout);

 private:
  int _verbosity;

};

#endif
