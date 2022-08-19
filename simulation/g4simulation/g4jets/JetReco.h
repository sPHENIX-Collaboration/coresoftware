#ifndef G4JET_JETRECO_H
#define G4JET_JETRECO_H

//===========================================================
/// \file JetReco.h
/// \brief simple jet reco using FastJet
/// \author Mike McCumber
//===========================================================

// PHENIX includes
#include <fun4all/SubsysReco.h>

// standard includes
#include <string>  // for string
#include <vector>

// forward declarations
class Jet;
class JetAlgo;
class JetInput;
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
  ~JetReco() override;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;

  void add_input(JetInput *input) { _inputs.push_back(input); }
  void add_algo(JetAlgo *algo, std::string output)
  {
    _algos.push_back(algo);
    _outputs.push_back(output);
  }

  void set_algo_node(const std::string &algonode) { _algonode = algonode; }
  void set_input_node(const std::string &inputnode) { _inputnode = inputnode; }

 private:
  int CreateNodes(PHCompositeNode *topNode);
  void FillJetNode(PHCompositeNode *topNode, int ialgo, std::vector<Jet *> jets);

  std::vector<JetInput *> _inputs;
  std::vector<JetAlgo *> _algos;
  std::string _algonode;
  std::string _inputnode;
  std::vector<std::string> _outputs;
};

#endif  // G4JET_JETRECO_H
