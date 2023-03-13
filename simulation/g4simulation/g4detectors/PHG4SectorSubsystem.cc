#include "PHG4SectorSubsystem.h"
#include "PHG4SectorDetector.h"
#include "PHG4SectorDisplayAction.h"
#include "PHG4SectorSteppingAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4Subsystem.h>       // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <sstream>

class PHG4Detector;

//_______________________________________________________________________
PHG4SectorSubsystem::PHG4SectorSubsystem(const std::string& name)
  : PHG4Subsystem(name)
{
}

//_______________________________________________________________________
PHG4SectorSubsystem::~PHG4SectorSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4SectorSubsystem::Init(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4SectorDisplayAction(Name());
  // create detector
  m_Detector = new PHG4SectorDetector(this, topNode, Name());
  m_Detector->geom = geom;
  m_Detector->SuperDetector(superdetector);
  m_Detector->OverlapCheck(CheckOverlap());

  if (geom.GetNumActiveLayers())
  {
    std::ostringstream nodename;
    if (superdetector != "NONE")
    {
      nodename << "G4HIT_" << superdetector;
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    // create hit list
    PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!block_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    // create stepping action
    m_SteppingAction = new PHG4SectorSteppingAction(m_Detector);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4SectorSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector*
PHG4SectorSubsystem::GetDetector() const
{
  return m_Detector;
}
