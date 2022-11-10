#include "BeamLineMagnetSubsystem.h"

#include "BeamLineMagnetDetector.h"
#include "BeamLineMagnetDisplayAction.h"
#include "BeamLineMagnetSteppingAction.h"
#include "PHG4DetectorSubsystem.h"  // for PHG4DetectorSubsystem

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <iostream>  // for operator<<, endl, basic_ostream
#include <set>       // for set

class PHG4Detector;

//_______________________________________________________________________
BeamLineMagnetSubsystem::BeamLineMagnetSubsystem(const std::string &na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
BeamLineMagnetSubsystem::~BeamLineMagnetSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int BeamLineMagnetSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  // create display settings before detector
  BeamLineMagnetDisplayAction *displayaction = new BeamLineMagnetDisplayAction(Name());
  m_DisplayAction = displayaction;
  /* create magnet */
  m_Detector = new BeamLineMagnetDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());

  if (GetParams()->get_int_param("active"))
  {
    std::set<std::string> nodes;
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dstNode;
    if (SuperDetector() != "NONE" || SuperDetector().empty())
    {
      // create super detector subnodes
      PHNodeIterator iter_dst(dstNode);
      DetNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(DetNode);
      }
    }
    // create hit output nodes
    std::string detector_suffix = SuperDetector();
    if (detector_suffix == "NONE")
    {
      detector_suffix = Name();
    }

    m_HitNodeName = "G4HIT_" + detector_suffix;
    nodes.insert(m_HitNodeName);
    m_AbsorberNodeName = "G4HIT_ABSORBER_" + detector_suffix;
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodes.insert(m_AbsorberNodeName);
    }
    for (const auto &node : nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node, "PHObject"));
      }
    }
    m_SteppingAction = new BeamLineMagnetSteppingAction(m_Detector, GetParams());
    m_SteppingAction->SetHitNodeName("G4HIT", m_HitNodeName);
    m_SteppingAction->SetHitNodeName("G4HIT_ABSORBER", m_AbsorberNodeName);
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new BeamLineMagnetSteppingAction(m_Detector, GetParams());
  }
  return 0;
}

//_______________________________________________________________________
int BeamLineMagnetSubsystem::process_event(PHCompositeNode *topNode)
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
PHG4Detector *BeamLineMagnetSubsystem::GetDetector() const
{
  return m_Detector;
}

void BeamLineMagnetSubsystem::SetDefaultParameters()
{
  set_default_string_param("magtype", "");

  set_default_double_param("field_x", 0.);
  set_default_double_param("field_y", 0.);
  set_default_double_param("field_z", 0.);
  set_default_double_param("fieldgradient", 0.);

  set_default_double_param("field_global_position_x", 0.);  // abs. position to world for field manager
  set_default_double_param("field_global_position_y", 0.);  // abs. position to world for field manager
  set_default_double_param("field_global_position_z", 0.);  // abs. position to world for field manager
  set_default_double_param("field_global_rot_x", 0.);       // abs. rotation to world for field manager
  set_default_double_param("field_global_rot_y", 0.);       // abs. rotation to world for field manager
  set_default_double_param("field_global_rot_z", 0.);       // abs. rotation to world for field manager

  set_default_double_param("length", 100);
  set_default_double_param("place_x", 0.);  // relative position to mother vol.
  set_default_double_param("place_y", 0.);  // relative position to mother vol.
  set_default_double_param("place_z", 0.);  // relative position to mother vol.
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("inner_radius", 4);
  set_default_double_param("outer_radius", 100);
  set_default_double_param("skin_thickness", 0.);  // Fe thickness before tracks are terminated
}

void BeamLineMagnetSubsystem::Print(const std::string & /*what*/) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  if (!BeginRunExecuted())
  {
    std::cout << "Need to execute BeginRun() before parameter printout is meaningful" << std::endl;
    std::cout << "To do so either run one or more events or on the command line execute: " << std::endl;
    std::cout << "Fun4AllServer *se = Fun4AllServer::instance();" << std::endl;
    std::cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << std::endl;
    std::cout << "g4->InitRun(se->topNode());" << std::endl;
    std::cout << "BeamLineMagnetSubsystem *cyl = (PHG4BeamLineMagnetSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << std::endl;
    std::cout << "cyl->Print()" << std::endl;
    return;
  }
  GetParams()->Print();
  return;
}
