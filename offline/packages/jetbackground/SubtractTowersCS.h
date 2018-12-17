#ifndef __SUBTRACTTOWERSCS_H__
#define __SUBTRACTTOWERSCS_H__

//===========================================================
/// \file SubtractTowersCS.h
/// \brief creates new UE-subtracted towers, using constituent subtraction
/// \author Dennis V. Perepelitsa
//===========================================================

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

// standard includes
#include <vector>

#include "calobase/RawTowerContainer.h"

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
  virtual ~SubtractTowersCS();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetFlowModulation(bool use_flow_modulation) { _use_flow_modulation = use_flow_modulation; }
  void SetAlpha(float alpha) { _alpha = alpha; }
  void SetDeltaRmax(float DeltaRmax) { _DeltaRmax = DeltaRmax; }

 private:
  int CreateNode(PHCompositeNode *topNode);

  bool _use_flow_modulation;

  float _alpha;
  float _DeltaRmax;
};

#endif  // __SUBTRACTTOWERSCS_H__
