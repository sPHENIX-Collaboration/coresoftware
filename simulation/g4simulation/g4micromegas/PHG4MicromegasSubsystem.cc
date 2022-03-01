/*!
 * \file PHG4MicromegasSubsystem.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "PHG4MicromegasSubsystem.h"

#include "PHG4MicromegasDetector.h"
#include "PHG4MicromegasDisplayAction.h"
#include "PHG4MicromegasSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

//_______________________________________________________________________
PHG4MicromegasSubsystem::PHG4MicromegasSubsystem(const std::string &name, int layerno)
  : PHG4DetectorSubsystem(name, layerno)
{
  // call base class method which will set up parameter infrastructure
  // and call our SetDefaultParameters() method
  InitializeParameters();

  SuperDetector(name);
}

//_______________________________________________________________________
PHG4MicromegasSubsystem::~PHG4MicromegasSubsystem()
{
  delete m_DisplayAction;
}


//_______________________________________________________________________
int PHG4MicromegasSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  PHNodeIterator dstIter(dstNode);

  if (GetParams()->get_int_param("active"))
  {
    auto detNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!detNode)
    {
      detNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(detNode);
    }

    const std::string g4hitnodename = "G4HIT_" + SuperDetector();
    auto g4_hits = findNode::getClass<PHG4HitContainer>(detNode, g4hitnodename);
    if (!g4_hits)
    {
      g4_hits = new PHG4HitContainer(g4hitnodename);
      detNode->addNode(new PHIODataNode<PHObject>(g4_hits, g4hitnodename, "PHObject"));
    }
  }

  // create detector
  m_DisplayAction = new PHG4MicromegasDisplayAction(Name());
  m_Detector = new PHG4MicromegasDetector(this, topNode, GetParams(), Name());
  m_Detector->set_first_layer( GetLayer() );

  m_Detector->Verbosity(Verbosity());
  m_Detector->SuperDetector(SuperDetector());
  m_Detector->OverlapCheck(CheckOverlap());
  
  // create stepping action if detector is active
  if (GetParams()->get_int_param("active"))
  { m_SteppingAction = new PHG4MicromegasSteppingAction(m_Detector, GetParams()); }
  return 0;
}

//_______________________________________________________________________
int PHG4MicromegasSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  { m_SteppingAction->SetInterfacePointers(topNode); }
  return 0;
}

//_______________________________________________________________________
void PHG4MicromegasSubsystem::Print(const std::string &what) const
{ if (m_Detector) m_Detector->Print(what); }

//_______________________________________________________________________
PHG4Detector* PHG4MicromegasSubsystem::GetDetector(void) const
{ return m_Detector; }

//_______________________________________________________________________
PHG4SteppingAction* PHG4MicromegasSubsystem::GetSteppingAction() const
{ return m_SteppingAction; }

//_______________________________________________________________________
void PHG4MicromegasSubsystem::SetDefaultParameters()
{}
