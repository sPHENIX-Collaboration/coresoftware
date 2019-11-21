#include "G4JLeicBeamLineMagnetSubsystem.h"
#include "G4JLeicBeamLineMagnetDetector.h"
#include "G4JLeicBeamLineMagnetDisplayAction.h"
#include "G4JLeicBeamLineMagnetSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>    // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>   // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>          // for PHIODataNode
#include <phool/PHNode.h>                // for PHNode
#include <phool/PHNodeIterator.h>        // for PHNodeIterator
#include <phool/PHObject.h>              // for PHObject
#include <phool/getClass.h>

#include <cmath>     // for NAN
#include <iostream>  // for operator<<, endl, basic_ostream

class PHCompositeNode;
class PHG4Detector;

using namespace std;

//_______________________________________________________________________
G4JLeicBeamLineMagnetSubsystem::G4JLeicBeamLineMagnetSubsystem(const std::string& na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
{
  InitializeParameters();
}

//_______________________________________________________________________
G4JLeicBeamLineMagnetSubsystem::~G4JLeicBeamLineMagnetSubsystem()
{
  delete  m_DisplayAction;
}

//_______________________________________________________________________
int G4JLeicBeamLineMagnetSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{

  // create display settings before detector
  G4JLeicBeamLineMagnetDisplayAction *displayaction = new G4JLeicBeamLineMagnetDisplayAction(Name(), GetParams());
  m_DisplayAction = displayaction;
  /* create magnet */
  m_Detector = new G4JLeicBeamLineMagnetDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));

    ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      // create super detector subnodes
      PHNodeIterator iter_dst(dstNode);
      PHCompositeNode *superSubNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(superSubNode);
      }
      dstNode = superSubNode;
      PHNodeIterator iter_run(runNode);
      superSubNode = dynamic_cast<PHCompositeNode *>(iter_run.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        runNode->addNode(superSubNode);
      }
      runNode = superSubNode;

      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    PHG4HitContainer *magnet_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!magnet_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(magnet_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    m_SteppingAction = new G4JLeicBeamLineMagnetSteppingAction(m_Detector, GetParams());
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new G4JLeicBeamLineMagnetSteppingAction(m_Detector, GetParams());
  }

  return 0;
}

//_______________________________________________________________________
int G4JLeicBeamLineMagnetSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

//_______________________________________________________________________
PHG4Detector* G4JLeicBeamLineMagnetSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void G4JLeicBeamLineMagnetSubsystem::SetDefaultParameters()
{
  set_default_string_param("magtype", "");

  set_default_double_param("field_x", 0.);
  set_default_double_param("field_y", 0.);
  set_default_double_param("field_z", 0.);
  set_default_double_param("fieldgradient", 0.);

  set_default_double_param("length", 100);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("inner_radius", 4);
  set_default_double_param("outer_radius", 100);
  set_default_string_param("material", "G4_Galactic");
}

void G4JLeicBeamLineMagnetSubsystem::Print(const string& what) const
{
  cout << Name() << " Parameters: " << endl;
  if (!BeginRunExecuted())
  {
    cout << "Need to execute BeginRun() before parameter printout is meaningful" << endl;
    cout << "To do so either run one or more events or on the command line execute: " << endl;
    cout << "Fun4AllServer *se = Fun4AllServer::instance();" << endl;
    cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << endl;
    cout << "g4->InitRun(se->topNode());" << endl;
    cout << "G4JLeicBeamLineMagnetSubsystem *cyl = (PHG4BeamLineMagnetSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << endl;
    cout << "cyl->Print()" << endl;
    return;
  }
  GetParams()->Print();
  return;
}
