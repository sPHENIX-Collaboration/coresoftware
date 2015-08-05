#ifndef __JETRECO_H__
#define __JETRECO_H__

//===========================================================
/// \file JetReco.h
/// \brief simple jet reco using FastJet
/// \author Mike McCumber
//===========================================================

#include "JetInput.h"
#include "JetAlgo.h"

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHTimeServer.h>

// standard includes
#include <vector>

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
  virtual ~JetReco();
		
  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int CreateNodes(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

  void add_input(JetInput* input) {_inputs.push_back(input);}
  void add_algo(JetAlgo* algo, std::string output) {
    _algos.push_back(algo);
    _outputs.push_back(output);
  }
  
 private:

  std::vector<JetInput*>  _inputs;
  std::vector<JetAlgo*>   _algos;
  std::vector<std::string> _outputs;
  
};

#endif // __JETRECO_H__
