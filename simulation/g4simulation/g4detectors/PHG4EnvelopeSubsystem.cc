#include "PHG4EnvelopeSubsystem.h"

#include "PHG4EnvelopeDetector.h"
#include "PHG4EnvelopeSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Subsystem.h>  // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <sstream>

using namespace std;

PHG4EnvelopeSubsystem::PHG4EnvelopeSubsystem(const std::string& name, const int /*lyr*/)
  : PHG4Subsystem(name)
  , detector_(nullptr)
  , steppingAction_(nullptr)
  , material("G4_PbWO4")
  ,  // default - lead tungstate crystal
  active(1)
  , detector_type(name)
{
}

int PHG4EnvelopeSubsystem::Init(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4EnvelopeDetector(this, topNode, Name());
  detector_->SetActive(active);
  detector_->OverlapCheck(CheckOverlap());

  if (active)
  {
    // create hit output node
    ostringstream nodename;
    nodename << "G4HIT_ENVELOPE_" << detector_type;

    PHG4HitContainer* crystal_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!crystal_hits)
    {
      crystal_hits = new PHG4HitContainer(nodename.str());
      PHIODataNode<PHObject>* hitNode = new PHIODataNode<PHObject>(crystal_hits, nodename.str(), "PHObject");
      dstNode->addNode(hitNode);
    }

    // create stepping action
    steppingAction_ = new PHG4EnvelopeSteppingAction(detector_);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4EnvelopeSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4EnvelopeSubsystem::GetDetector() const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4EnvelopeSubsystem::GetSteppingAction() const
{
  return steppingAction_;
}
