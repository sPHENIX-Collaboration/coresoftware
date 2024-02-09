#include "JetProbeMaker.h"

#include "Jetv2.h"
#include "JetContainer.h"
#include "JetContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/phool.h>  // for PHWHERE
#include <phool/getClass.h>

#include <TRandom3.h>
#include <fastjet/PseudoJet.hh>

int JetProbeMaker::process_event(PHCompositeNode* /*topNode*/) {
  // update the jet probe
  float phi = _rand->Uniform(_phi_min, _phi_max);
  float eta = _rand->Uniform(_eta_min, _eta_max);
  float pt = (_const_pt ? _pt : _rand->Uniform(_pt_min, _pt));

  fastjet::PseudoJet fjet {};
  fjet.reset_PtYPhiM(pt, eta, phi);

  Jetv2* jet = (Jetv2*) _jets->add_jet();
  jet->set_px(fjet.px());
  jet->set_py(fjet.py());
  jet->set_pz(fjet.pz());
  jet->set_e(fjet.e());
  jet->insert_comp(Jet::SRC::JET_PROBE, 0, false);

  return Fun4AllReturnCodes::EVENT_OK;
}

JetProbeMaker::JetProbeMaker(const std::string &name, const int _randseed) 
  : SubsysReco(name)
  , _rand { new TRandom3(_randseed)  }
{};

int JetProbeMaker::InitRun(PHCompositeNode *topNode) {
  // Create the Input node if required
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Input node if required
  PHCompositeNode *probeNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "JetProbeNode"));
  if (!probeNode)
  {
    probeNode = new PHCompositeNode("JetProbeNode");
    dstNode->addNode(probeNode);
  }

  _jets = findNode::getClass<JetContainer>(topNode, "JetProbeContainer");
  if (!_jets)
  {
    _jets = new JetContainerv1();
    PHIODataNode<PHObject> *JetContainerNode = new
      PHIODataNode<PHObject>(_jets, "JetProbeContainer", "PHObject");
    dstNode->addNode(JetContainerNode);
  }


  return Fun4AllReturnCodes::EVENT_OK;
}

JetProbeMaker::~JetProbeMaker () {
  delete _rand;
}
