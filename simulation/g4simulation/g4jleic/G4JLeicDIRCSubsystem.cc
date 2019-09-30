#include "G4JLeicDIRCSubsystem.h"
#include "G4JLeicDIRCDetector.h"
#include "G4JLeicDIRCSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>          // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                 // for PHIODataNode
#include <phool/PHNode.h>                       // for PHNode
#include <phool/PHNodeIterator.h>               // for PHNodeIterator
#include <phool/PHObject.h>                     // for PHObject
#include <phool/getClass.h>

#include <boost/foreach.hpp>

#include <set>                                  // for set
#include <sstream>

using namespace std;

//_______________________________________________________________________
G4JLeicDIRCSubsystem::G4JLeicDIRCSubsystem(const std::string &name)
  : PHG4DetectorSubsystem(name, 0)
  , detector_(nullptr)
  , steppingAction_(nullptr)
{
  InitializeParameters();
  Name(name);
  SuperDetector(name);
}

//_______________________________________________________________________
int G4JLeicDIRCSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new G4JLeicDIRCDetector(topNode, GetParams(), Name());
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());

  set<string> nodes;
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    nodes.insert(nodename.str());
    BOOST_FOREACH (string node, nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node.c_str());
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node.c_str(), "PHObject"));
      }
    }
    // create stepping action
    steppingAction_ = new G4JLeicDIRCSteppingAction(detector_, GetParams());
  }

  return 0;
}

//_______________________________________________________________________
int G4JLeicDIRCSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

void G4JLeicDIRCSubsystem::Print(const string &what) const
{
  //cout << "PSTOF Parameters: " << endl;
  if (detector_)
  {
    detector_->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *G4JLeicDIRCSubsystem::GetDetector(void) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction *G4JLeicDIRCSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
}

void G4JLeicDIRCSubsystem::SetDefaultParameters()
{
// all units are in cm
  set_default_double_param("PixelDy", 4./10.); // dy/10

  set_default_int_param( "layers", 6);
}
