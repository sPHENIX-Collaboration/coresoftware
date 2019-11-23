#include "G4JLeicVTXSubsystem.h"
#include "G4JLeicVTXDetector.h"
#include "G4JLeicVTXSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4detectors/PHG4DetectorGroupSubsystem.h>  // for PHG4DetectorGrou...

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <boost/foreach.hpp>

#include <set>  // for set
#include <sstream>

using namespace std;

//_______________________________________________________________________
G4JLeicVTXSubsystem::G4JLeicVTXSubsystem(const std::string &name)
  : PHG4DetectorGroupSubsystem(name)
  , detector_(nullptr)
  , steppingAction_(nullptr)
{
  InitializeParameters();
  Name(name);
  SuperDetector(name);
}

//_______________________________________________________________________
int G4JLeicVTXSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new G4JLeicVTXDetector(this, topNode, GetParamsContainer(), Name());
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());

  set<string> nodes;
  if (GetParamsContainer()->GetParameters(-1)->get_int_param("active"))
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
    steppingAction_ = new G4JLeicVTXSteppingAction(detector_, GetParamsContainer());
    steppingAction_->Init();
  }

  return 0;
}

//_______________________________________________________________________
int G4JLeicVTXSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

void G4JLeicVTXSubsystem::Print(const string &what) const
{
  //cout << "PSTOF Parameters: " << endl;
  PrintDefaultParams();
  PrintMacroParams();
  GetParamsContainer()->Print();
  if (detector_)
  {
    detector_->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *G4JLeicVTXSubsystem::GetDetector(void) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction *G4JLeicVTXSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
}

void G4JLeicVTXSubsystem::SetDefaultParameters()
{
  // all units are in cm
  set_default_double_param(0, "Dx", 0.005);
  set_default_double_param(0, "Dy", 2.);
  set_default_double_param(0, "Dz", 10.);
  set_default_double_param(0, "Rin", 3.5);
  set_default_double_param(0, "PixelDx", 1.);  // dz/10. = 1.
  set_default_double_param(0, "PixelDy", 2. / 50.);   // dy/50

  set_default_double_param(1, "Dx", 0.005);
  set_default_double_param(1, "Dy", 2.);
  set_default_double_param(1, "Dz", 11.);
  set_default_double_param(1, "Rin", 4.5);
  set_default_double_param(1, "PixelDx", 11. / 10.);  // dz/10
  set_default_double_param(1, "PixelDy", 2. / 50.);   // dy/50

  set_default_double_param(2, "Dx", 0.005);
  set_default_double_param(2, "Dy", 4.);
  set_default_double_param(2, "Dz", 18.);
  set_default_double_param(2, "Rin", 6.5);
  set_default_double_param(2, "PixelDx", 18. / 50.);  // dz/50
  set_default_double_param(2, "PixelDy", 4. / 10.);   // dy/10

  set_default_double_param(3, "Dx", 0.005);
  set_default_double_param(3, "Dy", 4.);
  set_default_double_param(3, "Dz", 24.);
  set_default_double_param(3, "Rin", 10.5);
  set_default_double_param(3, "PixelDx", 24. / 50.);  // dz/50
  set_default_double_param(3, "PixelDy", 4. / 10.);   // dy/10

  set_default_double_param(4, "Dx", 0.005);
  set_default_double_param(4, "Dy", 4.);
  set_default_double_param(4, "Dz", 36.);
  set_default_double_param(4, "Rin", 13.5);
  set_default_double_param(4, "PixelDx", 36. / 50.);  // dz/50
  set_default_double_param(4, "PixelDy", 4. / 10.);   // dy/10

  set_default_double_param(5, "Dx", 0.005);
  set_default_double_param(5, "Dy", 4.);
  set_default_double_param(5, "Dz", 48.);
  set_default_double_param(5, "Rin", 15.5);
  set_default_double_param(5, "PixelDx", 48. / 50.);  // dz/50
  set_default_double_param(5, "PixelDy", 4. / 10.);   // dy/10

  // geometry version number
  // we use negative numbers until the "official" version
  // when we build the detector
  // set_default_int_param(-1,"geometry_version",-1);
  set_default_int_param(-1, "layers", 6);
  set_default_int_param(-1, "rows", 56);
  set_default_double_param(-1, "shift_z", 40);
  set_default_double_param(-1, "xsize", 0.8);
  set_default_double_param(-1, "ysize", 6.);
  set_default_double_param(-1, "zsize", 5.);
  set_default_int_param(-1, "active", 1);
  set_default_int_param(-1, "absorberactive", 0);
}
