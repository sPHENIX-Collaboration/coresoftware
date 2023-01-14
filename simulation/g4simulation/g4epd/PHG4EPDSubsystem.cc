/* vim: set sw=2 ft=cpp: */

#include "PHG4EPDSubsystem.h"

#include "PHG4EPDDetector.h"
#include "PHG4EPDDisplayAction.h"
#include "PHG4EPDSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4detectors/PHG4DetectorSubsystem.h>  // for PHG4DetectorSubsystem

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <set>
#include <string>

PHG4EPDSubsystem::PHG4EPDSubsystem(std::string const& name)
  : PHG4DetectorSubsystem(name)
{
  InitializeParameters();
}

PHG4EPDSubsystem::~PHG4EPDSubsystem()
{
  delete m_DisplayAction;
}

int PHG4EPDSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  m_DisplayAction = new PHG4EPDDisplayAction(Name());

  m_Detector = new PHG4EPDDetector(this, topNode, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    std::set<std::string> nodes;
    PHNodeIterator iter(topNode);
    PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* DetNode = dstNode;
    if (SuperDetector() != "NONE" && !SuperDetector().empty())
    {
      PHNodeIterator iter_dst(dstNode);
      DetNode = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));

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

    for (const auto& nodename : nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4EPDSteppingAction(m_Detector, GetParams());
    m_SteppingAction->SetHitNodeName("G4HIT", m_HitNodeName);
    m_SteppingAction->SetHitNodeName("G4HIT_SUPPORT", m_SupportNodeName);
  }
  else if (GetParams()->get_int_param("blackhole"))
  {
    m_SteppingAction = new PHG4EPDSteppingAction(m_Detector, GetParams());
  }

  return 0;
}

int PHG4EPDSubsystem::process_event(PHCompositeNode* topNode)
{
  if (m_SteppingAction != nullptr)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

PHG4Detector* PHG4EPDSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4EPDSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_z", 316.);
}
