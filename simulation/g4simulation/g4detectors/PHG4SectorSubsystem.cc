#include "PHG4SectorSubsystem.h"
#include "PHG4SectorDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4SectorSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4SectorSubsystem::PHG4SectorSubsystem(const std::string &name) :
    PHG4Subsystem(name), detector_(0), steppingAction_(NULL), eventAction_(
        NULL), superdetector("NONE")
{

}

//_______________________________________________________________________
int
PHG4SectorSubsystem::Init(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4SectorDetector(topNode, Name());
  detector_->geom = geom;
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(CheckOverlap());

  if (geom.GetNumActiveLayers())
    {
      ostringstream nodename;
      if (superdetector != "NONE")
        {
          nodename << "G4HIT_" << superdetector;
        }
      else
        {
          nodename << "G4HIT_" << Name();
        }
      // create hit list
      PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(
          topNode, nodename.str().c_str());
      if (!block_hits)
        {

          dstNode->addNode(new PHIODataNode<PHObject>(block_hits =
						      new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject"));

        }
      // create stepping action
      steppingAction_ = new PHG4SectorSteppingAction(detector_);

      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4SectorSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector*
PHG4SectorSubsystem::GetDetector(void) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction*
PHG4SectorSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
}


