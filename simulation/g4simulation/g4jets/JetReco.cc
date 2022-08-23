
#include "JetReco.h"

#include "Jet.h"
#include "JetAlgo.h"
#include "JetInput.h"
#include "JetMap.h"
#include "JetMapv1.h"

// PHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// standard includes
#include <cstdlib>  // for exit
#include <iostream>
#include <memory>  // for allocator_traits<>::value_type
#include <vector>

JetReco::JetReco(const std::string &name)
  : SubsysReco(name)
{
}

JetReco::~JetReco()
{
  for (unsigned int i = 0; i < _inputs.size(); ++i)
  {
    delete _inputs[i];
  }
  _inputs.clear();
  for (unsigned int i = 0; i < _algos.size(); ++i)
  {
    delete _algos[i];
  }
  _algos.clear();
  _outputs.clear();
}

int JetReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "========================== JetReco::InitRun() =============================" << std::endl;
    std::cout << " Input Selections:" << std::endl;
    for (unsigned int i = 0; i < _inputs.size(); ++i) _inputs[i]->identify();
    std::cout << " Algorithms:" << std::endl;
    for (unsigned int i = 0; i < _algos.size(); ++i) _algos[i]->identify();
    std::cout << "===========================================================================" << std::endl;
  }

  return CreateNodes(topNode);
}

int JetReco::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 1) std::cout << "JetReco::process_event -- entered" << std::endl;

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  std::vector<Jet *> inputs;  // owns memory
  for (unsigned int iselect = 0; iselect < _inputs.size(); ++iselect)
  {
    std::vector<Jet *> parts = _inputs[iselect]->get_input(topNode);
    for (unsigned int ipart = 0; ipart < parts.size(); ++ipart)
    {
      inputs.push_back(parts[ipart]);
      inputs.back()->set_id(inputs.size() - 1);  // unique ids ensured
    }
  }

  //---------------------------
  // Run the jet reconstruction
  //---------------------------
  for (unsigned int ialgo = 0; ialgo < _algos.size(); ++ialgo)
  {
    std::vector<Jet *> jets = _algos[ialgo]->get_jets(inputs);  // owns memory

    // send the output somewhere on the DST
    FillJetNode(topNode, ialgo, jets);
  }

  // clean up input vector
  for (unsigned int i = 0; i < inputs.size(); ++i) delete inputs[i];
  inputs.clear();

  if (Verbosity() > 1) std::cout << "JetReco::process_event -- exited" << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the AntiKt node if required
  PHCompositeNode *AlgoNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _algonode));
  if (!AlgoNode)
  {
    AlgoNode = new PHCompositeNode(_algonode);
    dstNode->addNode(AlgoNode);
  }

  // Create the Input node if required
  PHCompositeNode *InputNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _inputnode));
  if (!InputNode)
  {
    InputNode = new PHCompositeNode(_inputnode);
    AlgoNode->addNode(InputNode);
  }

  for (unsigned i = 0; i < _outputs.size(); ++i)
  {
    JetMap *jets = findNode::getClass<JetMap>(topNode, _outputs[i]);
    if (!jets)
    {
      jets = new JetMapv1();
      PHIODataNode<PHObject> *JetMapNode = new PHIODataNode<PHObject>(jets, _outputs[i], "PHObject");
      InputNode->addNode(JetMapNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void JetReco::FillJetNode(PHCompositeNode *topNode, int ipos, std::vector<Jet *> jets)
{
  JetMap *jetmap = findNode::getClass<JetMap>(topNode, _outputs[ipos]);
  if (!jetmap)
  {
    std::cout << PHWHERE << " ERROR: Can't find JetMap: " << _outputs[ipos] << std::endl;
    exit(-1);
  }

  jetmap->set_algo(_algos[ipos]->get_algo());
  jetmap->set_par(_algos[ipos]->get_par());
  for (unsigned int i = 0; i < _inputs.size(); ++i)
  {
    jetmap->insert_src(_inputs[i]->get_src());
  }

  for (unsigned int i = 0; i < jets.size(); ++i)
  {
    jetmap->insert(jets[i]);  // map takes ownership, sets unique id
  }

  return;
}
