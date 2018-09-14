#ifndef __SUBTRACTTOWERS_H__
#define __SUBTRACTTOWERS_H__

//===========================================================
/// \file SubtractTowers.h
/// \brief creates new UE-subtracted towers
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

/// \class SubtractTowers
///
/// \brief creates new UE-subtracted towers
///
/// Using a previously determined background UE density, this module
/// constructs a new set of towers by subtracting the background from
/// existing raw towers
///
class SubtractTowers : public SubsysReco
{
 public:
  SubtractTowers(const std::string &name = "SubtractTowers");
  virtual ~SubtractTowers();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void SetFlowModulation( bool use_flow_modulation ) { _use_flow_modulation = use_flow_modulation; }

 private:
  int CreateNode(PHCompositeNode *topNode);

  int _background_iteration;

  bool _use_flow_modulation;
};

#endif  // __SUBTRACTTOWERS_H__
