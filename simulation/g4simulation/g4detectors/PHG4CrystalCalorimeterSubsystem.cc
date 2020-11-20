#include "PHG4CrystalCalorimeterSubsystem.h"
#include "PHG4CrystalCalorimeterDetector.h"
#include "PHG4CrystalCalorimeterDisplayAction.h"
#include "PHG4CrystalCalorimeterSteppingAction.h"
#include "PHG4ProjCrystalCalorimeterDetector.h"

#include <g4main/PHG4DisplayAction.h>              // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>             // for PHG4SteppingAction
#include <g4main/PHG4Subsystem.h>                  // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                    // for PHIODataNode
#include <phool/PHNode.h>                          // for PHNode
#include <phool/PHNodeIterator.h>                  // for PHNodeIterator
#include <phool/PHObject.h>                        // for PHObject
#include <phool/getClass.h>

#include <iostream>                                // for operator<<, ostrin...
#include <sstream>

class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4CrystalCalorimeterSubsystem::PHG4CrystalCalorimeterSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
  , active(1)
  , detector_type(name)
  , mappingfile_("")
  , mappingfile_4x4_construct_("")
  , projective_(false)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4CrystalCalorimeterSubsystem::~PHG4CrystalCalorimeterSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4CrystalCalorimeterSubsystem::InitSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4CrystalCalorimeterDisplayAction(Name());
  // create detector
  if (projective_)
  {
    if (Verbosity() > 1)
    {
      cout << "PHG4CrystalCalorimeterSubsystem::InitRun - use PHG4ProjCrystalCalorimeterDetector" << endl;
    }
    m_Detector = new PHG4ProjCrystalCalorimeterDetector(this, topNode, Name());
    m_Detector->SetTowerMappingFile(mappingfile_);
    m_Detector->SetSupermoduleGeometry(mappingfile_4x4_construct_);
  }
  else
  {
    if (Verbosity() > 1)
    {
      cout << "PHG4CrystalCalorimeterSubsystem::InitRun - use PHG4CrystalCalorimeterDetector" << endl;
    }
    m_Detector = new PHG4CrystalCalorimeterDetector(this, topNode, Name());
    m_Detector->SetTowerMappingFile(mappingfile_);
  }

  m_Detector->SetActive(active);
  m_Detector->SetAbsorberActive(active);
  m_Detector->OverlapCheck(CheckOverlap());

  if (active)
  {
    // create hit output node
    ostringstream nodename;
    nodename << "G4HIT_" << detector_type;

    PHG4HitContainer* crystal_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
    if (!crystal_hits)
    {
      crystal_hits = new PHG4HitContainer(nodename.str());
      PHIODataNode<PHObject>* hitNode = new PHIODataNode<PHObject>(crystal_hits, nodename.str().c_str(), "PHObject");
      dstNode->addNode(hitNode);
    }

    ostringstream absnodename;
    absnodename << "G4HIT_ABSORBER_" << detector_type;

    PHG4HitContainer* absorber_hits = findNode::getClass<PHG4HitContainer>(topNode, absnodename.str().c_str());
    if (!absorber_hits)
    {
      absorber_hits = new PHG4HitContainer(absnodename.str());
      PHIODataNode<PHObject>* abshitNode = new PHIODataNode<PHObject>(absorber_hits, absnodename.str().c_str(), "PHObject");
      dstNode->addNode(abshitNode);
    }

    // create stepping action
    m_SteppingAction = new PHG4CrystalCalorimeterSteppingAction(m_Detector);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4CrystalCalorimeterSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4CrystalCalorimeterSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4CrystalCalorimeterSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", -108.);
  return;
}
