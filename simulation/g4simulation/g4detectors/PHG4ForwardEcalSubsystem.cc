#include "PHG4ForwardEcalSubsystem.h"
#include "PHG4EICForwardEcalDetector.h"
#include "PHG4ForwardEcalDetector.h"
#include "PHG4ForwardEcalDisplayAction.h"
#include "PHG4ForwardEcalSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <TSystem.h>

#include <cstdlib>                         // for getenv
#include <set>  // for set
#include <sstream>

class PHG4Detector;

//_______________________________________________________________________
PHG4ForwardEcalSubsystem::PHG4ForwardEcalSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4ForwardEcalSubsystem::~PHG4ForwardEcalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4ForwardEcalSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4ForwardEcalDisplayAction(Name());
  // create detector
  if (m_EICDetectorFlag)
  {
    m_Detector = new PHG4EICForwardEcalDetector(this, topNode, GetParams(), Name());
  }
  else
  {
    m_Detector = new PHG4ForwardEcalDetector(this, topNode, GetParams(), Name());
  }

  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());
  std::set<std::string> nodes;
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* DetNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    std::string nodename;
    if (SuperDetector() != "NONE")
    {
      nodename = "G4HIT_" + SuperDetector();
    }
    else
    {
      nodename = "G4HIT_" + Name();
    }
    nodes.insert(nodename);

    if (GetParams()->get_int_param("absorberactive"))
    {
      if (SuperDetector() != "NONE")
      {
        nodename = "G4HIT_ABSORBER_" + SuperDetector();
      }
      else
      {
        nodename = "G4HIT_ABSORBER_" + Name();
      }
      nodes.insert(nodename);
    }

    for (auto thisnode : nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, thisnode);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(thisnode);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, thisnode, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4ForwardEcalSteppingAction(m_Detector, GetParams());
  }

  return 0;
}

//_______________________________________________________________________
int PHG4ForwardEcalSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4ForwardEcalSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4ForwardEcalSubsystem::SetDefaultParameters()
{
  std::ostringstream mappingfilename;
  const char* calibroot = getenv("CALIBRATIONROOT");
  if (calibroot)
  {
    mappingfilename << calibroot;
  }
  else
  {
    std::cout << "no CALIBRATIONROOT environment variable" << std::endl;
    gSystem->Exit(1);
  }
  mappingfilename << "/ForwardEcal/mapping/towerMap_FEMC_fsPHENIX_v004.txt";
  set_default_string_param("mapping_file", mappingfilename.str());
  set_default_string_param("mapping_file_md5", PHG4Utils::md5sum(mappingfilename.str()));
  return;
}

void PHG4ForwardEcalSubsystem::SetTowerMappingFile(const std::string& filename)
{
  set_string_param("mapping_file", filename);
  set_string_param("mapping_file_md5", PHG4Utils::md5sum(get_string_param("mapping_file")));
}
