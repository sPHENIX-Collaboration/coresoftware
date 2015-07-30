#ifndef __JETRECO_H__
#define __JETRECO_H__

//===========================================================
/// \file JetReco.h
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

/// \class JetReco
///
/// \brief jet reco with user def inputs and algos
///
/// This module can be used to reconstruct truth jets
/// and will get me started on filling some jet nodes and getting
/// source material for jet evaluation
///
class JetReco : public SubsysReco
{

 public:
 
  JetReco(const std::string &name = "JetReco");
  virtual ~JetReco() {}
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 private:

};

#endif // __JETRECO_H__
