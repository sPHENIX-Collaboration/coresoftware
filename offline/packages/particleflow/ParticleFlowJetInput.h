#ifndef PARTICLEFLOW_PARTICLEFLOWJETINPUT_H
#define PARTICLEFLOW_PARTICLEFLOWJETINPUT_H

//===========================================================
/// \file ParticleFlowJetInput.h
/// \brief Connective tissue between jet reco and PFlow elements
/// \author Dennis V. Perepelitsa
//===========================================================

#include <g4jets/JetInput.h>

// finally system includes
#include <iostream>    // for cout, ostream


// forward declarations
class PHCompositeNode;
class Jet;

class ParticleFlowJetInput : public JetInput
{
 public:
  ParticleFlowJetInput() = default;
  ~ParticleFlowJetInput() override {}

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;
  void identify(std::ostream& os = std::cout) override;

 private:
  int _verbosity = 0;

};

#endif
