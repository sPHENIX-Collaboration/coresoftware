#include "PHG4FCalSubsystem.h"

#include "PHG4FCalDetector.h"
#include "PHG4FCalSteppingAction.h"

#include <g4main/PHG4HitContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4UserSteppingAction.hh>  // for G4UserSteppingAction

class PHG4Detector;

//_______________________________________________________________________
PHG4FCalSubsystem::PHG4FCalSubsystem(const char* name)
  : PHG4Subsystem(name)
  , detector_(0)
{
}

//_______________________________________________________________________
int PHG4FCalSubsystem::Init(PHCompositeNode* topNode)
{
  // create hit list
  PHG4HitContainer* fcal_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_FCAL");
  if (!fcal_hits)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    dstNode->addNode(new PHIODataNode<PHObject>(fcal_hits = new PHG4HitContainer("G4HIT_FCAL"), "G4HIT_FCAL", "PHObject"));
  }

  // create detector
  detector_ = new PHG4FCalDetector(this, topNode, Name());
  detector_->Verbosity(Verbosity());

  // create stepping action

  return 9;
}

//_______________________________________________________________________
int PHG4FCalSubsystem::process_event(PHCompositeNode* topNode)
{
  if (PHG4FCalSteppingAction* p = dynamic_cast<PHG4FCalSteppingAction*>(detector_->GetSteppingAction()))
    p->SetInterfacePointers(topNode);

  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4FCalSubsystem::GetDetector(void) const
{
  return detector_;
}
