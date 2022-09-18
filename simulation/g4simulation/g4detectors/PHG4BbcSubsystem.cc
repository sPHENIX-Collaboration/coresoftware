#include "PHG4BbcSubsystem.h"

#include "PHG4BbcDetector.h"
#include "PHG4BbcDisplayAction.h"
#include "PHG4BbcSteppingAction.h"

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
#include <phool/phool.h>

#include <iostream>  // for operator<<, basic_ostream, endl
#include <set>       // for set
#include <sstream>

//_______________________________________________________________________
PHG4BbcSubsystem::PHG4BbcSubsystem(const std::string &name)
  : PHG4DetectorSubsystem(name)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4BbcSubsystem::~PHG4BbcSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4BbcSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    std::cout << PHWHERE << " Initializing BBC Subsystem" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4BbcDisplayAction(Name());
  // create detector
  m_Detector = new PHG4BbcDetector(this, topNode, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    std::set<std::string> nodes;
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dstNode;
    if (SuperDetector() != "NONE" && !SuperDetector().empty())
    {
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
    m_SupportNodeName = "G4HIT_SUPPORT_" + detector_suffix;
    if (GetParams()->get_int_param("supportactive"))
    {
      nodes.insert(m_SupportNodeName);
    }
    for (const auto &nodename : nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4BbcSteppingAction(m_Detector, GetParams());
    m_SteppingAction->SetHitNodeName("G4HIT", m_HitNodeName);
    m_SteppingAction->SetHitNodeName("G4HIT_SUPPORT", m_SupportNodeName);
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new PHG4BbcSteppingAction(m_Detector, GetParams());
  }

  return 0;
}

//_______________________________________________________________________
int PHG4BbcSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4BbcSubsystem::Print(const std::string &what) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  GetParams()->Print();
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  if (m_SteppingAction)
  {
    m_SteppingAction->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4BbcSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4BbcSubsystem::SetDefaultParameters()
{
  // geometry version number
  // we use negative numbers until the "official" version
  // when we build the detector
  // set_default_int_param("geometry_version",-1);
}
