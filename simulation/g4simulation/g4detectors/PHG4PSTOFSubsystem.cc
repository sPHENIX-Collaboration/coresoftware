// $Id$

/*!
 * \file ${file_name}
 * \brief
 * \author Mickey Chiu <chiu@bnl.gov>
 * \version $Revision$
 * \date $Date$
 */

#include "PHG4PSTOFSubsystem.h"
#include "PHG4PSTOFDetector.h"
#include "PHG4PSTOFSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>  // for PHG4SteppingAction

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <boost/foreach.hpp>

#include <set>  // for set
#include <sstream>

//_______________________________________________________________________
PHG4PSTOFSubsystem::PHG4PSTOFSubsystem(const std::string &name)
  : PHG4DetectorGroupSubsystem(name)
{
  InitializeParameters();
  Name(name);
  SuperDetector(name);
}

//_______________________________________________________________________
int PHG4PSTOFSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4PSTOFDetector(this, topNode, GetParamsContainer(), Name());
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());

  if (GetParamsContainer()->GetParameters(-1)->get_int_param("active"))
  {
    std::set<std::string> nodes;
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }
    std::ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    nodes.insert(nodename.str());
    BOOST_FOREACH (std::string node, nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node);
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node, "PHObject"));
      }
    }
    // create stepping action
    steppingAction_ = new PHG4PSTOFSteppingAction(detector_, GetParamsContainer());
    steppingAction_->Init();
  }

  return 0;
}

//_______________________________________________________________________
int PHG4PSTOFSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4PSTOFSubsystem::Print(const std::string &what) const
{
  //std::cout << "PSTOF Parameters: " << std::endl;
  PrintDefaultParams();
  PrintMacroParams();
  GetParamsContainer()->Print();
  if (detector_)
  {
    detector_->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4PSTOFSubsystem::GetDetector() const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction *PHG4PSTOFSubsystem::GetSteppingAction() const
{
  return steppingAction_;
}

void PHG4PSTOFSubsystem::SetDefaultParameters()
{
  set_default_double_param(0, "z_mod_0", -109.3);
  set_default_double_param(1, "z_mod_0", -96.66);
  set_default_double_param(2, "z_mod_0", -84.42);
  set_default_double_param(3, "z_mod_0", -72.55);
  set_default_double_param(4, "z_mod_0", -61.07);
  set_default_double_param(5, "z_mod_0", -49.97);
  set_default_double_param(6, "z_mod_0", -39.25);
  set_default_double_param(7, "z_mod_0", -28.72);
  set_default_double_param(8, "z_mod_0", -18.76);
  set_default_double_param(9, "z_mod_0", -9.191);
  set_default_double_param(10, "z_mod_0", 0);
  set_default_double_param(11, "z_mod_0", 9.191);
  set_default_double_param(12, "z_mod_0", 18.76);
  set_default_double_param(13, "z_mod_0", 28.72);
  set_default_double_param(14, "z_mod_0", 39.25);
  set_default_double_param(15, "z_mod_0", 49.97);
  set_default_double_param(16, "z_mod_0", 61.07);
  set_default_double_param(17, "z_mod_0", 72.55);
  set_default_double_param(18, "z_mod_0", 84.42);
  set_default_double_param(19, "z_mod_0", 96.66);
  set_default_double_param(20, "z_mod_0", 109.3);

  set_default_double_param(0, "z_mod_1", -107.2);
  set_default_double_param(1, "z_mod_1", -94.66);
  set_default_double_param(2, "z_mod_1", -82.52);
  set_default_double_param(3, "z_mod_1", -70.75);
  set_default_double_param(4, "z_mod_1", -59.37);
  set_default_double_param(5, "z_mod_1", -48.47);
  set_default_double_param(6, "z_mod_1", -37.85);
  set_default_double_param(7, "z_mod_1", -27.72);
  set_default_double_param(8, "z_mod_1", -18.76);
  set_default_double_param(9, "z_mod_1", -9.191);
  set_default_double_param(10, "z_mod_1", 0);
  set_default_double_param(11, "z_mod_1", 9.191);
  set_default_double_param(12, "z_mod_1", 18.76);
  set_default_double_param(13, "z_mod_1", 27.72);
  set_default_double_param(14, "z_mod_1", 37.85);
  set_default_double_param(15, "z_mod_1", 48.47);
  set_default_double_param(16, "z_mod_1", 59.37);
  set_default_double_param(17, "z_mod_1", 70.75);
  set_default_double_param(18, "z_mod_1", 82.52);
  set_default_double_param(19, "z_mod_1", 94.66);
  set_default_double_param(20, "z_mod_1", 107.2);

  set_default_double_param(0, "r_mod_0", 85.6);
  set_default_double_param(1, "r_mod_0", 85.6);
  set_default_double_param(2, "r_mod_0", 85.6);
  set_default_double_param(3, "r_mod_0", 85.6);
  set_default_double_param(4, "r_mod_0", 86);
  set_default_double_param(5, "r_mod_0", 86.5);
  set_default_double_param(6, "r_mod_0", 86.5);
  set_default_double_param(7, "r_mod_0", 86.5);
  set_default_double_param(8, "r_mod_0", 85.5);
  set_default_double_param(9, "r_mod_0", 83.6);
  set_default_double_param(10, "r_mod_0", 87.5);
  set_default_double_param(11, "r_mod_0", 83.6);
  set_default_double_param(12, "r_mod_0", 85.5);
  set_default_double_param(13, "r_mod_0", 86.5);
  set_default_double_param(14, "r_mod_0", 86.5);
  set_default_double_param(15, "r_mod_0", 86.5);
  set_default_double_param(16, "r_mod_0", 86);
  set_default_double_param(17, "r_mod_0", 85.6);
  set_default_double_param(18, "r_mod_0", 85.6);
  set_default_double_param(19, "r_mod_0", 85.6);
  set_default_double_param(20, "r_mod_0", 85.6);

  set_default_double_param(0, "r_mod_1", 85.3);
  set_default_double_param(1, "r_mod_1", 85.2);
  set_default_double_param(2, "r_mod_1", 84.9);
  set_default_double_param(3, "r_mod_1", 84.8);
  set_default_double_param(4, "r_mod_1", 85.1);
  set_default_double_param(5, "r_mod_1", 85);
  set_default_double_param(6, "r_mod_1", 85);
  set_default_double_param(7, "r_mod_1", 84.8);
  set_default_double_param(8, "r_mod_1", 83.8);
  set_default_double_param(9, "r_mod_1", 81.9);
  set_default_double_param(10, "r_mod_1", 85.8);
  set_default_double_param(11, "r_mod_1", 81.9);
  set_default_double_param(12, "r_mod_1", 83.8);
  set_default_double_param(13, "r_mod_1", 84.8);
  set_default_double_param(14, "r_mod_1", 85);
  set_default_double_param(15, "r_mod_1", 85);
  set_default_double_param(16, "r_mod_1", 85.1);
  set_default_double_param(17, "r_mod_1", 84.8);
  set_default_double_param(18, "r_mod_1", 84.9);
  set_default_double_param(19, "r_mod_1", 85.2);
  set_default_double_param(20, "r_mod_1", 85.3);

  // geometry version number
  // we use negative numbers until the "official" version
  // when we build the detector
  // set_default_int_param(-1,"geometry_version",-1);
  set_default_int_param(-1, "modules", 21);
  set_default_int_param(-1, "rows", 56);
  set_default_double_param(-1, "xsize", 0.8);
  set_default_double_param(-1, "ysize", 6.);
  set_default_double_param(-1, "zsize", 5.);
  set_default_int_param(-1, "active", 1);
  set_default_int_param(-1, "absorberactive", 0);
}
