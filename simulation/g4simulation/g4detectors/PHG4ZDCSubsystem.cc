#include "PHG4ZDCSubsystem.h"
#include "PHG4ZDCDetector.h"
#include "PHG4ZDCDisplayAction.h"
#include "PHG4ZDCSteppingAction.h"

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

#include <cstdlib>                         // for getenv
#include <set>  // for set
#include <fstream>
#include <iostream>
#include <sstream>

class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4ZDCSubsystem::PHG4ZDCSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
 
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4ZDCSubsystem::~PHG4ZDCSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4ZDCSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector
  m_DisplayAction = new PHG4ZDCDisplayAction(Name());
  // create detector
 
  m_Detector = new PHG4ZDCDetector(this, topNode, GetParams(), Name());

  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  m_Detector->Verbosity(Verbosity());
  set<string> nodes;
  // GetParams()->set_int_param("active",1);
  // std::cout<<"active?"<<GetParams()->get_int_param("active")<<std::endl;
  
  if (GetParams()->get_int_param("active"))
   {
   
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode* DetNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    ostringstream nodename;
    
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    nodes.insert(nodename.str());

    if (GetParams()->get_int_param("absorberactive"))
    {
      nodename.str("");
      if (SuperDetector() != "NONE")
      {
        nodename << "G4HIT_ABSORBER_" << SuperDetector();
      }
      else
      {
        nodename << "G4HIT_ABSORBER_" << Name();
      }
      nodes.insert(nodename.str());
    }
    for (auto nodename : nodes)

    //    BOOST_FOREACH (string node, nodes)
    {
      PHG4HitContainer* g4_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(nodename);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, nodename, "PHObject"));
      }
    }
    // create stepping action
    m_SteppingAction = new PHG4ZDCSteppingAction(m_Detector, GetParams());
  }

  return 0;
}

//_______________________________________________________________________
int PHG4ZDCSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4ZDCSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4ZDCSubsystem::SetDefaultParameters()
{
  
  return;
}

