#include "PHG4SectorSubsystem.h"
#include "PHG4SectorDetector.h"
#include "PHG4SectorDisplayAction.h"
#include "PHG4SectorSteppingAction.h"

#include <g4main/PHG4HitContainer.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4SectorSubsystem::PHG4SectorSubsystem(const std::string& name)
  : PHG4Subsystem(name)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
  , superdetector("NONE")
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
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4SectorDisplayAction(Name());
  // create detector
  m_Detector = new PHG4SectorDetector(this, topNode, Name());
  m_Detector->geom = geom;
  m_Detector->SuperDetector(superdetector);
  m_Detector->OverlapCheck(CheckOverlap());

  if (geom.GetNumActiveLayers())
  {
    ostringstream nodename;
    if (superdetector != "NONE")
    {
      nodename << "G4HIT_" << superdetector;
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    // create hit list
    PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(
        topNode, nodename.str().c_str());
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
PHG4SectorSubsystem::GetDetector(void) const
{
  return m_Detector;
}
