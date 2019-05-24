#include "PHG4ForwardHcalSubsystem.h"
#include "PHG4ForwardHcalDetector.h"
#include "PHG4ForwardHcalDisplayAction.h"
#include "PHG4ForwardHcalSteppingAction.h"

#include <g4main/PHG4HitContainer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4ForwardHcalSubsystem::PHG4ForwardHcalSubsystem(const std::string& name, const int lyr)
  : PHG4Subsystem(name)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
  , active(1)
  , absorber_active(0)
  , blackhole(0)
  , detector_type(name)
  , mappingfile_("")
{
}

//_______________________________________________________________________
PHG4ForwardHcalSubsystem::~PHG4ForwardHcalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4ForwardHcalSubsystem::Init(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4ForwardHcalDisplayAction(Name());
  // create detector
  m_Detector = new PHG4ForwardHcalDetector(this, topNode, Name());
  m_Detector->SetActive(active);
  m_Detector->SetAbsorberActive(absorber_active);
  m_Detector->BlackHole(blackhole);
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());
  m_Detector->SetTowerMappingFile(mappingfile_);

  if (active)
  {
    set<string> nodes;

    // create hit output node
    ostringstream nodename;
    nodename << "G4HIT_" << detector_type;

    PHG4HitContainer* scintillator_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
    if (!scintillator_hits)
    {
      scintillator_hits = new PHG4HitContainer(nodename.str());
      PHIODataNode<PHObject>* hitNode = new PHIODataNode<PHObject>(scintillator_hits, nodename.str().c_str(), "PHObject");
      dstNode->addNode(hitNode);
      nodes.insert(nodename.str());
    }

    ostringstream absnodename;
    absnodename << "G4HIT_ABSORBER_" << detector_type;

    PHG4HitContainer* absorber_hits = findNode::getClass<PHG4HitContainer>(topNode, absnodename.str().c_str());
    if (!absorber_hits)
    {
      absorber_hits = new PHG4HitContainer(absnodename.str());
      PHIODataNode<PHObject>* abshitNode = new PHIODataNode<PHObject>(absorber_hits, absnodename.str().c_str(), "PHObject");
      dstNode->addNode(abshitNode);
      nodes.insert(nodename.str());
    }

    // create stepping action
    m_SteppingAction = new PHG4ForwardHcalSteppingAction(m_Detector);
  }

  return 0;
}

//_______________________________________________________________________
int PHG4ForwardHcalSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4ForwardHcalSubsystem::GetDetector() const
{
  return m_Detector;
}
