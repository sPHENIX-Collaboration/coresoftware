#ifndef JETBACKGROUND_SUBTRACTTOWERSCS_H
#define JETBACKGROUND_SUBTRACTTOWERSCS_H

//===========================================================
/// \file SubtractTowersCS.h
/// \brief creates new UE-subtracted towers, using constituent subtraction
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <string> 

// forward declarations
class PHCompositeNode;

/// \class SubtractTowersCS
///
/// \brief creates new UE-subtracted towers
///
/// Using a previously determined background UE density, this module
/// constructs a new set of towers by subtracting the background from
/// existing raw towers. CS parameters are configurable
///
class SubtractTowersCS : public SubsysReco
{
 public:
  SubtractTowersCS(const std::string &name = "SubtractTowersCS");
  ~SubtractTowersCS() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void SetFlowModulation(bool use_flow_modulation) { _use_flow_modulation = use_flow_modulation; }
  void SetAlpha(float alpha) { _alpha = alpha; }
  void SetDeltaRmax(float DeltaRmax) { _DeltaRmax = DeltaRmax; }

 private:
  int CreateNode(PHCompositeNode *topNode);

  bool _use_flow_modulation;

  float _alpha;
  float _DeltaRmax;
};

#endif
