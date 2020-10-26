// $Id$

/*!
 * \file ${file_name}
 * \brief
 * \author Mickey Chiu <chiu@bnl.gov>
 * \version $Revision$
 * \date $Date$
 */

#include "PHG4BbcSubsystem.h"
#include "PHG4BbcDetector.h"
#include "PHG4BbcSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>          // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                 // for PHIODataNode
#include <phool/PHNode.h>                       // for PHNode
#include <phool/PHNodeIterator.h>               // for PHNodeIterator
#include <phool/PHObject.h>                     // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <boost/foreach.hpp>

#include <set>                                  // for set
#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4BbcSubsystem::PHG4BbcSubsystem(const std::string &name)
  : PHG4DetectorSubsystem(name)
  , m_detector(nullptr)
  , m_steppingAction(nullptr)
{
  InitializeParameters();
  Name(name);
  SuperDetector(name);
}

//_______________________________________________________________________
int PHG4BbcSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  cout << PHWHERE << endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  m_detector = new PHG4BbcDetector(this, topNode, GetParams(), Name());
  m_detector->SuperDetector(SuperDetector());
  m_detector->OverlapCheck(CheckOverlap());

  set<string> nodes;
  //if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
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
    BOOST_FOREACH (string node, nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node.c_str());
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node.c_str(), "PHObject"));
      }
    }
    // create stepping action
    m_steppingAction = new PHG4BbcSteppingAction(m_detector, GetParams());
    m_steppingAction->Init();
  }

  return 0;
}

//_______________________________________________________________________
int PHG4BbcSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_steppingAction)
  {
    m_steppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4BbcSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  GetParams()->Print();
  if (m_detector)
  {
    m_detector->Print(what);
  }
  if (m_steppingAction)
  {
    m_steppingAction->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4BbcSubsystem::GetDetector(void) const
{
  return m_detector;
}

//_______________________________________________________________________
PHG4SteppingAction *PHG4BbcSubsystem::GetSteppingAction(void) const
{
  return m_steppingAction;
}

void PHG4BbcSubsystem::SetDefaultParameters()
{
  set_default_double_param("z", 250.);

  // geometry version number
  // we use negative numbers until the "official" version
  // when we build the detector
  // set_default_int_param("geometry_version",-1);
  //set_default_int_param("active", 1);
  //set_default_int_param("absorberactive", 0);
}
