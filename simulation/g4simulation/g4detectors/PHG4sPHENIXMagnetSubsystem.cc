#include "PHG4sPHENIXMagnetSubsystem.h"

#include "PHG4HcalDefs.h"
#include "PHG4sPHENIXMagnetDetector.h"
#include "PHG4sPHENIXMagnetDisplayAction.h"
#include "PHG4sPHENIXMagnetSteppingAction.h"

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

#include <cmath>     // for NAN
#include <iostream>  // for operator<<, basic_ostream
#include <set>       // for set
#include <sstream>

class PHG4Detector;

//_______________________________________________________________________
PHG4sPHENIXMagnetSubsystem::PHG4sPHENIXMagnetSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4sPHENIXMagnetSubsystem::~PHG4sPHENIXMagnetSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4sPHENIXMagnetSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4sPHENIXMagnetDisplayAction(Name());

  // create detector
  m_Detector = new PHG4sPHENIXMagnetDetector(this, topNode, GetParams(), Name(), GetLayer());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    std::set<std::string> nodes;
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    std::ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    nodes.insert(nodename.str());

    for (auto &node : nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4sPHENIXMagnetSteppingAction(m_Detector, GetParams());
    m_SteppingAction->InitWithNode(topNode);
  }
  else
  {
    if (GetParams()->get_int_param("blackhole"))
    {
      m_SteppingAction = new PHG4sPHENIXMagnetSteppingAction(m_Detector, GetParams());
      m_SteppingAction->InitWithNode(topNode);
    }
  }

  return 0;
}

//_______________________________________________________________________
int PHG4sPHENIXMagnetSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4sPHENIXMagnetSubsystem::Print(const std::string &what) const
{
  std::cout << "sPHENIX Magnet Parameters: " << std::endl;
  GetParams()->Print();
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4sPHENIXMagnetSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4sPHENIXMagnetSubsystem::SetDefaultParameters()
{
  return;
}
