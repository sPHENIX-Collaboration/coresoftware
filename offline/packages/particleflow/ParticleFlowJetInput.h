#ifndef PARTICLEFLOW_PARTICLEFLOWJETINPUT_H
#define PARTICLEFLOW_PARTICLEFLOWJETINPUT_H

//===========================================================
/// \file ParticleFlowJetInput.h
/// \brief Connective tissue between jet reco and PFlow elements
/// \author Dennis V. Perepelitsa
//===========================================================

// finally system includes
#include <iostream>    // for cout, ostream

#include <g4jets/Jet.h>
#include <g4jets/JetInput.h>

// forward declarations
class PHCompositeNode;

class ParticleFlowJetInput : public JetInput
{
 public:
  ParticleFlowJetInput();
  ~ParticleFlowJetInput() override {}

  std::vector<Jet*> get_input(PHCompositeNode* topNode) override;
  void identify(std::ostream& os = std::cout) override;

 private:
  int _verbosity;

};

#endif
