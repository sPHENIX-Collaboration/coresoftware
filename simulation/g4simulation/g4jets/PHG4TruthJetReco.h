#ifndef __PHG4TRUTHJETRECO_H__
#define __PHG4TRUTHJETRECO_H__

//===========================================================
/// \file PHG4TruthJetReco.h
/// \brief simple jet reco using FastJet
/// \author Mike McCumber
//===========================================================

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimeServer.h>

// standard includes

// forward declarations
class PHCompositeNode;

/// \class PHG4TruthJetReco
///
/// \brief simple jet reco using FastJet
///
/// This module can be used to reconstruct truth jets
/// and will get me started on filling some jet nodes and getting
/// source material for jet evaluation
///
class PHG4TruthJetReco : public SubsysReco
{

 public:
 
  PHG4TruthJetReco(const std::string &name = "PHG4TruthJetReco");
  virtual ~PHG4TruthJetReco() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

};

#endif // __PHG4TRUTHJETRECO_H__
