#include "PHG4SiliconTrackerSubsystem.h"
#include "PHG4Parameters.h"
#include "PHG4ParametersContainer.h"
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
  , superdetector("NONE")
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

  int active = 0;
  int absorberactive = 0;
  int blackhole = 0;
  for (set<int>::const_iterator parcontaineriter = GetDetIds().first; parcontaineriter != GetDetIds().second; ++parcontaineriter)
  {
    if (active || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("active"))
    {
      active = 1;
    }
    if (absorberactive || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("absorberactive"))
    {
      absorberactive = 1;
    }
    if (blackhole || GetParamsContainer()->GetParameters(*parcontaineriter)->get_int_param("blackhole"))
    {
      blackhole = 1;
    }
  }
  if (active)
  {
    cout << "detector: " << detector_type << endl;
    cout << "superdetector: " << superdetector << endl;
    std::string nodename = (superdetector != "NONE") ? boost::str(boost::format("G4HIT_%s") % superdetector) : boost::str(boost::format("G4HIT_%s") % detector_type);

    // create hit list
    PHG4HitContainer *hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
    if (!hitcontainer)
      dstNode->addNode(new PHIODataNode<PHObject>(hitcontainer = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));

    if (absorberactive)
    {
      //      nodename = (superdetector != "NONE") ? boost::str(boost::format("G4HIT_ABSORBER_%s") % superdetector) : boost::str(boost::format("G4HIT_ABSORBER_%s") % detector_type);

      hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, nodename.c_str());
      if (!hitcontainer)
      {
        dstNode->addNode(new PHIODataNode<PHObject>(hitcontainer = new PHG4HitContainer(nodename), nodename.c_str(), "PHObject"));
      }
    }

    // create stepping action
    steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParamsContainer());
  }
  else
  {
    if (blackhole)
    {
      steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_, GetParamsContainer());
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
// all values in cm!
  set_default_int_param(0,"nladder",20);
  set_default_int_param(1,"nladder",26);
  set_default_int_param(2,"nladder",32);
  set_default_int_param(3,"nladder",38);

  set_default_double_param(0, "radius", 6.);
  set_default_double_param(1, "radius", 8.);
  set_default_double_param(2, "radius", 10.);
  set_default_double_param(3, "radius", 12.);

  set_default_double_param(0, "strip_y",0.0078);
  set_default_double_param(1, "strip_y",0.0086);
  set_default_double_param(2, "strip_y",0.0086);
  set_default_double_param(3, "strip_y",0.0086);

  set_default_double_param(0, "strip_z_0",1.8);
  set_default_double_param(1, "strip_z_0",1.6);
  set_default_double_param(2, "strip_z_0",1.6);
  set_default_double_param(3, "strip_z_0",1.6);

  set_default_double_param(0, "strip_z_1",1.8);
  set_default_double_param(1, "strip_z_1",2.0);
  set_default_double_param(2, "strip_z_1",2.0);
  set_default_double_param(3, "strip_z_1",2.0);

  std::pair<std::set<int>::const_iterator, std::set<int>::const_iterator> begin_end = GetDetIds();
  for (set<int>::const_iterator it = begin_end.first; it != begin_end.second; ++it)
  {
    set_default_double_param(*it,"strip_x",0.02);
    set_int_param(*it,"active",1);
  }

  return;
}

void PHG4SiliconTrackerSubsystem::Print(const string &what) const
{
  PrintDefaultParams();
  PrintMacroParams();
  cout << "Print Resulting Parameters:" << endl;
  GetParamsContainer()->Print();
}
