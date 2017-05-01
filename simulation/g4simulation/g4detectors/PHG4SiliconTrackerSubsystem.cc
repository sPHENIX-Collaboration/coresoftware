#include "PHG4SiliconTrackerSubsystem.h"
#include "PHG4Parameters.h"
#include "PHG4SiliconTrackerDetector.h"
#include "PHG4SiliconTrackerSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

#include <boost/format.hpp>

using namespace std;

//_______________________________________________________________________
PHG4SiliconTrackerSubsystem::PHG4SiliconTrackerSubsystem(const std::string &detectorname, const vpair &layerconfig)
  : PHG4DetectorGroupSubsystem(detectorname)
  , detector_(0)
  , steppingAction_(nullptr)
  , layerconfig_(layerconfig)
  , detector_type(detectorname)
  , superdetector(detectorname)
{
  for (vector<pair<int, int>>::const_iterator piter = layerconfig.begin(); piter != layerconfig.end(); ++piter)
  {
    AddDetId((*piter).second);
  }
  InitializeParameters();
  // put the layer into the name so we get unique names
  // for multiple layers
  Name(detectorname);
}

//_______________________________________________________________________
int PHG4SiliconTrackerSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  if (verbosity > 0)
    std::cout << "PHG4SiliconTrackerSubsystem::Init started" << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4SiliconTrackerDetector(topNode, GetParamsContainer(), Name(), layerconfig_);
  detector_->SuperDetector(superdetector);
  detector_->Detector(detector_type);
  detector_->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    const int sphxlayermin = layerconfig_.front().first;
    const int sphxlayermax = layerconfig_.back().first;

    std::string nodename = (superdetector != "NONE") ? boost::str(boost::format("G4HIT_%s") % superdetector) : boost::str(boost::format("G4HIT_%s_%d_%d") % detector_type % sphxlayermin % sphxlayermax);

    // create hit list
    PHG4HitContainer *block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
    if (!block_hits)
      dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));

//    if (GetParams()->get_int_param("absorberactive"))
    if (0)
    {
      nodename = (superdetector != "NONE") ? boost::str(boost::format("G4HIT_ABSORBER_%s") % superdetector) : boost::str(boost::format("G4HIT_ABSORBER_%s_%d_%d") % detector_type % sphxlayermin % sphxlayermax);

      block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
      if (!block_hits)
        dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));

    }

    // create stepping action
    steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParams());
  }
  else
  {
    if (GetParams()->get_int_param("blackhole"))
    {
      steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParams());
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4SiliconTrackerSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
    steppingAction_->SetInterfacePointers(topNode);

  return 0;
}

//_______________________________________________________________________
PHG4Detector *PHG4SiliconTrackerSubsystem::GetDetector(void) const
{
  return detector_;
}

void PHG4SiliconTrackerSubsystem::SetDefaultParameters()
{
  set_default_double_param(0,"Radius",6.);
  set_default_double_param(1,"Radius",8.);
  set_default_double_param(2,"Radius",10.);
  set_default_double_param(3,"Radius",12.);
  return;
}

void PHG4SiliconTrackerSubsystem::Print(const string &what) const
{
  PrintDefaultParams();
}
