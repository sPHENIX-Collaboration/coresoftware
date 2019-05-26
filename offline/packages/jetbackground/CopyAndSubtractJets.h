#ifndef JETBACKGROUND_COPYANDSUBTRACTJETS_H
#define JETBACKGROUND_COPYANDSUBTRACTJETS_H

//===========================================================
/// \file CopyAndSubtractJets.h
/// \brief Creates subtracted copy of a jet collection
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string>

// forward declarations
class PHCompositeNode;

/// \class CopyAndSubtractJets
///
/// \brief Creates subtractd copy of a jet collection
///
/// Makes a copy of a jet collection with a new name and then updates
/// the kinematics of the jets in that collection based on a given UE
/// background (intended use is to create the set of jets used as
/// seeds in the second part of UE determination procedure)
///
class CopyAndSubtractJets : public SubsysReco
{
 public:
  CopyAndSubtractJets(const std::string &name = "CopyAndSubtractJets");
  virtual ~CopyAndSubtractJets() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetFlowModulation(bool use_flow_modulation) { _use_flow_modulation = use_flow_modulation; }

 private:
  int CreateNode(PHCompositeNode *topNode);

  bool _use_flow_modulation;
};

#endif
