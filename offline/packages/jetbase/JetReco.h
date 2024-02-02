#ifndef JETBASE_JETRECO_H
#define JETBASE_JETRECO_H

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
   enum TRANSITION
   {
    JET_CONTAINER,
    JET_MAP,
    BOTH,
    PRETEND_BOTH // do just JET_CONTAINER, but still append _JC to the name
   };

  JetReco(const std::string &name = "JetReco", TRANSITION _which_fill=JET_CONTAINER); // , bool fill_JetContainer=false);
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
  /* void set_fill_JetContainer(bool b) { _fill_JetContainer = b; } */

  JetAlgo* get_algo(unsigned int which_algo=0);

 private:
  int CreateNodes(PHCompositeNode *topNode);
  void FillJetNode(PHCompositeNode *topNode, int ialgo, std::vector<Jet *> jets);
  void FillJetContainer(PHCompositeNode *topNode, int ialgo, std::vector<Jet*>& jets);

  std::vector<JetInput *> _inputs;
  std::vector<JetAlgo *> _algos;
  std::string _algonode;
  std::string _inputnode;
  std::vector<std::string> _outputs;

  // transition functions, while moving from JetMap to JetContainer.
  // May be removed after transition is made, depending on state of 
  // functions
  std::string JC_name (std::string name) { 
    if (which_fill == TRANSITION::BOTH || which_fill==TRANSITION::PRETEND_BOTH) return name+"_JC";
    else return name;
  }
  TRANSITION which_fill;// fill both container and map 
  bool use_jetcon;
  bool use_jetmap;
};

#endif  // JETBASE_JETRECO_H
