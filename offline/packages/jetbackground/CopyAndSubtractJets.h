#ifndef __COPYANDSUBTRACTJETS_H__
#define __COPYANDSUBTRACTJETS_H__

//===========================================================
/// \file CopyAndSubtractJets.h
/// \brief Creates subtracted copy of a jet collection
/// \author Dennis V. Perepelitsa
//===========================================================

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <phool/PHTimeServer.h>

// standard includes
#include <vector>

#include <calobase/RawTowerContainer.h>

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
  virtual ~CopyAndSubtractJets();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:
  int CreateNode(PHCompositeNode *topNode);

};

#endif  // __COPYANDSUBTRACTJETS_H__
