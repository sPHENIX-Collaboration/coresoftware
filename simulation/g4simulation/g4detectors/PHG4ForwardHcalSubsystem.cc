#include "PHG4ForwardHcalSubsystem.h"

#include "PHG4ForwardHcalDetector.h"
#include "PHG4ForwardHcalDisplayAction.h"
#include "PHG4ForwardHcalSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4DisplayAction.h>       // for PHG4DisplayAction
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>      // for PHG4SteppingAction
#include <g4main/PHG4Subsystem.h>           // for PHG4Subsystem
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>             // for PHIODataNode
#include <phool/PHNode.h>                   // for PHNode
#include <phool/PHNodeIterator.h>           // for PHNodeIterator
#include <phool/PHObject.h>                 // for PHObject
#include <phool/getClass.h>

#include <TSystem.h>

#include <set>                              // for set
#include <sstream>

class PHG4Detector;

//_______________________________________________________________________
PHG4ForwardHcalSubsystem::PHG4ForwardHcalSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4ForwardHcalSubsystem::~PHG4ForwardHcalSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4ForwardHcalSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4ForwardHcalDisplayAction(Name());
  // create detector
  m_Detector = new PHG4ForwardHcalDetector(this, topNode, GetParams(), Name());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());

  std::set<std::string> nodes;
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* DetNode = dstNode;
    if (SuperDetector() != "NONE")
    {
      DetNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
      if (!DetNode)
      {
	DetNode = new PHCompositeNode(SuperDetector());
	dstNode->addNode(DetNode);
      }
    }
    // create hit output node
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
    for (auto nodename : nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4ForwardHcalSteppingAction(m_Detector,GetParams());
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

void PHG4ForwardHcalSubsystem::SetDefaultParameters()
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

  mappingfilename << "/ForwardHcal/mapping/towerMap_FHCAL_v005.txt";
  set_default_string_param("mapping_file", mappingfilename.str());
  set_default_string_param("mapping_file_md5", PHG4Utils::md5sum(mappingfilename.str()));
  return;
}

void PHG4ForwardHcalSubsystem::SetTowerMappingFile(const std::string& filename)
{
  set_string_param("mapping_file", filename);
  set_string_param("mapping_file_md5", PHG4Utils::md5sum(get_string_param("mapping_file")));
}
