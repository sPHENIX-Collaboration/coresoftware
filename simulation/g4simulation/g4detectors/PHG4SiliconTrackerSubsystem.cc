#include "PHG4SiliconTrackerSubsystem.h"
#include "PHG4SiliconTrackerDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4SiliconTrackerSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

#include <boost/format.hpp>

//_______________________________________________________________________
PHG4SiliconTrackerSubsystem::PHG4SiliconTrackerSubsystem(const std::string &detectorname, const vpair &layerconfig):
    PHG4Subsystem(detectorname),
    detector_(0),
    steppingAction_(NULL),
    eventAction_(NULL),
    active(0),
    absorberactive(0),
    layerconfig_(layerconfig),
    detector_type(detectorname),
    superdetector(detectorname)
{
  // put the layer into the name so we get unique names
  // for multiple layers
  Name(detectorname);
  for (int i = 0; i < 3; i++)
    dimension[i] = 100.0 * cm;
}

//_______________________________________________________________________
int PHG4SiliconTrackerSubsystem::Init( PHCompositeNode* topNode )
{
  if(verbosity>0)
    std::cout << "PHG4SiliconTrackerSubsystem::Init started" << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4SiliconTrackerDetector(topNode, Name(), layerconfig_);
  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->SuperDetector(superdetector);
  detector_->Detector(detector_type);
  detector_->OverlapCheck(overlapcheck);

  if (active)
    {
      const int sphxlayermin = layerconfig_.front().first;
      const int sphxlayermax = layerconfig_.back().first;

      std::string nodename = (superdetector != "NONE") ? boost::str(boost::format("G4HIT_%s") %superdetector) : boost::str(boost::format("G4HIT_%s_%d_%d") %detector_type %sphxlayermin %sphxlayermax);

      // create hit list
      PHG4HitContainer* block_hits =  findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
      if ( !block_hits )
        dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));

      PHG4EventActionClearZeroEdep *eventaction = new PHG4EventActionClearZeroEdep(topNode, nodename);
      if (absorberactive)
        {
          nodename = (superdetector != "NONE") ? boost::str(boost::format("G4HIT_ABSORBER_%s") %superdetector) : boost::str(boost::format("G4HIT_ABSORBER_%s_%d_%d") %detector_type %sphxlayermin %sphxlayermax);

          block_hits =  findNode::getClass<PHG4HitContainer>(topNode , nodename.c_str());
          if (!block_hits)
            dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));

          eventaction->AddNode(nodename);
        }

      eventAction_ = dynamic_cast<PHG4EventAction *> (eventaction);

      // create stepping action
      steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_);
    }

  return 0;
}

//_______________________________________________________________________
int PHG4SiliconTrackerSubsystem::process_event(PHCompositeNode * topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
    steppingAction_->SetInterfacePointers(topNode);

  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4SiliconTrackerSubsystem::GetDetector(void) const
  {
    return detector_;
  }

//_______________________________________________________________________
PHG4SteppingAction* PHG4SiliconTrackerSubsystem::GetSteppingAction(void) const
  {
    return steppingAction_;
  }

