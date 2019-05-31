// $Id: $

/*!
 * \file PHG4GDMLSubsystem.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHG4GDMLSubsystem.h"

#include "PHG4GDMLDetector.h"

#include <phparameter/PHParameters.h>

#include <iostream>                    // for operator<<, basic_ostream, endl

class PHCompositeNode;
class PHG4Detector;

using namespace std;

PHG4GDMLSubsystem::PHG4GDMLSubsystem(const std::string &name)
  : PHG4DetectorSubsystem(name, 0)
  , m_Detector(nullptr)
{
  InitializeParameters();
}

PHG4GDMLSubsystem::~PHG4GDMLSubsystem()
{
}

//_______________________________________________________________________
int PHG4GDMLSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
//  PHNodeIterator iter(topNode);
//  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  m_Detector = new PHG4GDMLDetector(topNode, Name(), GetParams());
  m_Detector->OverlapCheck(CheckOverlap());

  //  set<string> nodes;
  //  if (GetParams()->get_int_param("active"))
  //  {
  //    PHNodeIterator dstIter(dstNode);
  //    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
  //    if (!DetNode)
  //    {
  //      DetNode = new PHCompositeNode(SuperDetector());
  //      dstNode->addNode(DetNode);
  //    }
  //
  //    ostringstream nodename;
  //    if (SuperDetector() != "NONE")
  //    {
  //      nodename << "G4HIT_" << SuperDetector();
  //    }
  //    else
  //    {
  //      nodename << "G4HIT_" << Name();
  //    }
  //    nodes.insert(nodename.str());
  //    if (GetParams()->get_int_param("absorberactive"))
  //    {
  //      nodename.str("");
  //      if (SuperDetector() != "NONE")
  //      {
  //        nodename << "G4HIT_ABSORBER_" << SuperDetector();
  //      }
  //      else
  //      {
  //        nodename << "G4HIT_ABSORBER_" << Name();
  //      }
  //      nodes.insert(nodename.str());
  //    }
  //    BOOST_FOREACH (string node, nodes)
  //    {
  //      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node.c_str());
  //      if (!g4_hits)
  //      {
  //        g4_hits = new PHG4HitContainer(node);
  //        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node.c_str(), "PHObject"));
  //      }
  //    }
  //
  //    // create stepping action
  //    m_SteppingAction = new PHG4InnerHcalSteppingAction(m_Detector, GetParams());
  //  }
  //  else
  //  {
  //    // if this is a black hole it does not have to be active
  //    if (GetParams()->get_int_param("blackhole"))
  //    {
  //      m_SteppingAction = new PHG4InnerHcalSteppingAction(m_Detector, GetParams());
  //    }
  //  }
  return 0;
}

//_______________________________________________________________________
int PHG4GDMLSubsystem::process_event(PHCompositeNode *topNode)
{
  return 0;
}

void PHG4GDMLSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  GetParams()->Print();
  if (m_Detector)
  {
    m_Detector->Print(what);
  }

  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4GDMLSubsystem::GetDetector(void) const
{
  return m_Detector;
}

void PHG4GDMLSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);

  set_default_string_param("GDMPath", "DefaultParameters-InvadPath");
  set_default_string_param("TopVolName", "DefaultParameters-InvadVol");
}
