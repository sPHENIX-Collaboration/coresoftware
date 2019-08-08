// This is a single layer subsystem
// If you want multiple layers, then add more of these
// The layer number only affects the name

#include "PHG4OuterTrackerSubsystem.h"

#include "PHG4OuterTrackerDefs.h"
#include "PHG4OuterTrackerDetector.h"
#include "PHG4OuterTrackerDisplayAction.h"
#include "PHG4OuterTrackerSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>


#include <g4detectors/PHG4DetectorGroupSubsystem.h>  // for PHG4DetectorGrou...

#include <g4main/PHG4DisplayAction.h>                // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>

#include <phool/PHIODataNode.h>                      // for PHIODataNode
#include <phool/PHNode.h>                            // for PHNode
#include <phool/PHNodeIterator.h>                    // for PHNodeIterator
#include <phool/PHObject.h>                          // for PHObject
#include <phool/phool.h>                             // for PHWHERE
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <cmath>                                    // for NAN
#include <iostream>                                  // for operator<<, basi...
#include <set>                                       // for _Rb_tree_const_i...
#include <sstream>
#include <utility>                                   // for pair

class PHG4Detector;
class PHG4SteppingAction;

using namespace std;

//_______________________________________________________________________
PHG4OuterTrackerSubsystem::PHG4OuterTrackerSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorGroupSubsystem(name, lyr)
  , m_Detector(nullptr)
  , steppingAction_(nullptr)
  , m_DisplayAction(nullptr)
  , layer(lyr)
  , detector_type(name)
{

  InitializeParameters();

  // put the layers into name so we get unique names
  // for multiple layers
  Name(name);
  SuperDetector(name);
}

//_______________________________________________________________________
PHG4OuterTrackerSubsystem::~PHG4OuterTrackerSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4OuterTrackerSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  if (Verbosity() > 0)
    cout << "PHG4OuterTrackerSubsystem::Init started" << endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector (detector adds its volumes to it)
  m_DisplayAction = new PHG4OuterTrackerDisplayAction(Name());
  // create detector
  // These values are set from the calling macro using the setters defined in the .h file
  if (Verbosity())
  {
    cout << "    create OuterTracker detector with layer number " << layer  << endl;
  }
  m_Detector = new PHG4OuterTrackerDetector(this, layer, topNode, GetParamsContainer(), Name());
  m_Detector->Verbosity(Verbosity());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->Detector(detector_type);
  m_Detector->OverlapCheck(CheckOverlap());
  if (Verbosity())
  {
    cout << "    ------ created detector " << Name() << endl;
    GetParamsContainer()->Print();
  }
  const PHParameters *params = GetParamsContainer()->GetParameters(layer);

  if (params->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* detNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!detNode)
    {
      detNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(detNode);
    }
    ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << detector_type;
    }
    // create hit list
    PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!block_hits)
    {
      detNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    if (Verbosity())
      cout << PHWHERE << "creating hits node " << nodename.str() << endl;

    if (params->get_int_param("absorberactive"))
    {
      nodename.str("");
      if (SuperDetector() != "NONE")
      {
        nodename << "G4HIT_ABSORBER_" << SuperDetector();
      }
      else
      {
        nodename << "G4HIT_ABSORBER_" << detector_type;
      }
      block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
      if (!block_hits)
      {
        detNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
      }
      if (Verbosity())
        cout << PHWHERE << "creating hits node " << nodename.str() << endl;
    }

    // create stepping action
    steppingAction_ = new PHG4OuterTrackerSteppingAction(m_Detector, layer);
  }
  else
  {
   if (params->get_int_param("blackhole"))
    {
      steppingAction_ = new PHG4OuterTrackerSteppingAction(m_Detector, layer);
    }
  }
  return 0;
}

//_______________________________________________________________________
int PHG4OuterTrackerSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4OuterTrackerSubsystem::GetDetector(void) const
{
  return m_Detector;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4OuterTrackerSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
}

void PHG4OuterTrackerSubsystem::SetDefaultParameters()
{

  set_default_int_param(layer, "active", 1);  //non-automatic initialization in PHG4DetectorGroupSubsystem
  set_default_int_param(layer, "absorberactive", 0);  //non-automatic initialization in PHG4DetectorGroupSubsystem
  set_default_int_param(layer, "blackhole", 0);  //non-automatic initialization in PHG4DetectorGroupSubsystem
  set_default_int_param(layer, "layer", 0);
  set_default_double_param(layer, "ot_inner_radius", 80.0);  // cm
  set_default_double_param(layer, "ot_outer_radius", 80.01);
  set_default_double_param(layer, "ot_length", 120.0);
  set_default_int_param(layer, "ot_nseg_phi", 5026);
  set_default_int_param(layer, "ot_nseg_z", 1200);

}
