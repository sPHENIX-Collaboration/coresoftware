
#include "JetReco.h"

#include "Jet.h"
#include "JetAlgo.h"
#include "JetContainer.h"
#include "JetContainerv1.h"
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
#include <fstream>

JetReco::JetReco(const std::string &name, TRANSITION _which)
  : SubsysReco(name)
  , which_fill{_which}
  , use_jetcon{_which == TRANSITION::JET_CONTAINER || _which == TRANSITION::BOTH || _which == TRANSITION::PRETEND_BOTH}
  , use_jetmap{_which == TRANSITION::JET_MAP || _which == TRANSITION::BOTH}
{
}
/* JetReco::JetReco(const std::string &name, TRANSITION _which) */
/*   : SubsysReco(name) */
/*   , which_fill{_which} */
/*   , use_jetcon{false} */
/*   , use_jetmap{true} */
/* { */
/* } */


JetReco::~JetReco()
{
  for (auto &_input : _inputs)
  {
    delete _input;
  }
  _inputs.clear();
  for (auto &_algo : _algos)
  {
    delete _algo;
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
    for (auto &_input : _inputs) _input->identify();
    std::cout << " Algorithms:" << std::endl;
    for (auto &_algo : _algos) _algo->identify();
    std::cout << "===========================================================================" << std::endl;
  }

  return CreateNodes(topNode);
}

int JetReco::process_event(PHCompositeNode *topNode)
{

  if (Verbosity() > 1) std::cout << "JetReco::process_event -- entered" << std::endl;

  //------------------------------------------------------------------
  // This will also need to go into TClonesArrays in a future revision
  // Get Objects off of the Node Tree
  //------------------------------------------------------------------

  std::vector<Jet *> inputs;  // owns memory
  for (auto &_input : _inputs)
  {
    std::vector<Jet *> parts = _input->get_input(topNode);
    for (auto &part : parts)
    {
      inputs.push_back(part);
      inputs.back()->set_id(inputs.size() - 1);  // unique ids ensured
    }
  }

  //---------------------------
  // Run the jet reconstruction
  //---------------------------
  for (unsigned int ialgo = 0; ialgo < _algos.size(); ++ialgo)
  {
    // send the output somewhere on the DST
    /* if (_fill_JetContainer) { */
    if (use_jetcon)
    {
        if (Verbosity() > 5) std::cout << " Verbosity>5:: filling JetContainter for " << JC_name(_outputs[ialgo]) << std::endl;
        FillJetContainer(topNode, ialgo, inputs);
    }
    if (use_jetmap)
    {
      if (Verbosity() > 5) std::cout << " Verbosity>5:: filling jetnode for " << _outputs[ialgo] << std::endl;
      std::vector<Jet *> jets = _algos[ialgo]->get_jets(inputs);  // owns memory
      FillJetNode(topNode, ialgo, jets);
    }

    if (false ) { // These printouts were used in comparing the two map methods
                  // It can be removed later
      if (use_jetmap) {
        std::fstream fout;
        fout.open("fixme_jetmap", std::fstream::app);
        fout << " Printing jet results " << std::endl;
        JetMap *jetmap = findNode::getClass<JetMap>(topNode, _outputs[ialgo]);
        int ifixme =0;
        for (auto _jet = jetmap->begin(); _jet != jetmap->end(); ++_jet) {
          auto jet = _jet->second;
          fout << Form(" jet[%2i] ncon:pt:eta:phi [%6i,%6.3f,%6.3f,%6.3f]",
              ++ifixme, (int)jet->size_comp(), jet->get_pt(), jet->get_eta(), jet->get_phi()) << std::endl;
          std::vector<std::pair<int,int>> vconst;
          /*legacy*/ for (auto _comp = jet->begin_comp(); _comp != jet->end_comp(); ++_comp) {
            vconst.push_back({(int)_comp->first,(int)_comp->second});
          }
          std::sort(vconst.begin(), vconst.end(), [](std::pair<int,int> a, std::pair<int,int> b) 
              { 
              if      (a.first == b.first) return a.second < b.second;
              else return a.first < b.first;
              });
          fout << " constituents: ";
          int iconst = 0;
          for (const auto& p : vconst) {
            fout << Form("  c(%3i) %3i->%3i", iconst++, p.first, p.second) << std::endl;
          }
          fout << std::endl;
        }
      }
      if (use_jetcon) {
        std::fstream fout;
        fout.open("fixme_jetcont", std::fstream::app);
        fout << " Printing jet results" << std::endl;
        JetContainer *jet_cont = findNode::getClass<JetContainer>(topNode,JC_name(_outputs[ialgo]));
        int ifixme =0;

        /* for (auto jet = jet_cont->begin(); jet != jet_cont->end(); ++jet) { */
        for (auto jet : *jet_cont) {
          fout << Form(" jet[%2i] ncon:pt:eta:phi [%6i,%6.3f,%6.3f,%6.3f]",
              ++ifixme, (int)jet->size_comp(), jet->get_pt(), jet->get_eta(), jet->get_phi()) << std::endl;

          fout << " constituents: ";
          int iconst = 0;
          for (const auto& p : jet->get_comp_vec()) {
            fout << Form("  c(%3i) %3i->%3i", iconst++, p.first, p.second) << std::endl;
          }

          fout << std::endl;
        }
      }
    }

  }

  // clean up input vector
  // <- another place where TClonesArray's would make this more efficient
  for (auto &input : inputs) delete input;
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

  for (auto &_output : _outputs)
  {
    if (use_jetcon)
    {
      JetContainer *jetconn = findNode::getClass<JetContainer>(topNode, JC_name(_output));
      if (!jetconn)
      {
        jetconn = new JetContainerv1();
        PHIODataNode<PHObject> *JetContainerNode = new PHIODataNode<PHObject>(jetconn, JC_name(_output), "PHObject");
        InputNode->addNode(JetContainerNode);
      }
    }
    if (use_jetmap)
    {
      JetMap *jets = findNode::getClass<JetMap>(topNode, _output);
      if (!jets)
      {
        jets = new JetMapv1();
        PHIODataNode<PHObject> *JetMapNode = new PHIODataNode<PHObject>(jets, _output, "PHObject");
        InputNode->addNode(JetMapNode);
      }
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
  for (auto &_input : _inputs)
  {
    jetmap->insert_src(_input->get_src());
  }

  for (auto &jet : jets)
  {
    jetmap->insert(jet);  // map takes ownership, sets unique id
  }
  
  return;
}

void JetReco::FillJetContainer(PHCompositeNode *topNode, int ipos, std::vector<Jet *> &inputs)
{
  JetContainer *jetconn = findNode::getClass<JetContainer>(topNode, JC_name(_outputs[ipos]));
  if (!jetconn)
  {
    std::cout << PHWHERE << " ERROR: Can't find JetContainer: " << _outputs[ipos] << std::endl;
    exit(-1);
  }
  _algos[ipos]->cluster_and_fill(inputs, jetconn);  // fills the jet container with clustered jets
  for (auto &_input : _inputs)
  {
    jetconn->insert_src(_input->get_src());
  }

  if (Verbosity() > 7)
  {
    std::cout << " Verbosity()>7:: jets in container " << _outputs[ipos] << std::endl;
    jetconn->print_jets();
  }

  return;
}

JetAlgo *JetReco::get_algo(unsigned int which_algo)
{
  if (_algos.size() == 0)
  {
    std::cout << PHWHERE << std::endl
              << " JetReco has only " << _algos.size() << " JetAlgos; cannot get the one indexed " << which_algo << std::endl;
    assert(which_algo < _algos.size());
    return nullptr;
  }
  return _algos[which_algo];
}
