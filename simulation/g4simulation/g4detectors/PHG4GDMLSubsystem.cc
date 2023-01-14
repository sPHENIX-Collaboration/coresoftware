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

#include <iostream>  // for operator<<, basic_ostream, endl

class PHCompositeNode;
class PHG4Detector;

PHG4GDMLSubsystem::PHG4GDMLSubsystem(const std::string &name)
  : PHG4DetectorSubsystem(name, 0)
{
  InitializeParameters();
}

PHG4GDMLSubsystem::~PHG4GDMLSubsystem() = default;

//_______________________________________________________________________
int PHG4GDMLSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  //  PHNodeIterator iter(topNode);
  //  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  m_Detector = new PHG4GDMLDetector(this, topNode, Name(), GetParams());
  m_Detector->OverlapCheck(CheckOverlap());

  return 0;
}

//_______________________________________________________________________
int PHG4GDMLSubsystem::process_event(PHCompositeNode * /*topNode*/)
{
  return 0;
}

void PHG4GDMLSubsystem::Print(const std::string &what) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  GetParams()->Print();
  if (m_Detector)
  {
    m_Detector->Print(what);
  }

  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4GDMLSubsystem::GetDetector() const
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

  set_default_int_param("skip_DST_geometry_export", 0.);
}
