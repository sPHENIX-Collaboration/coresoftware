#include "PHG4BeamlineMagnetSubsystem.h"
#include "PHG4BeamlineMagnetDetector.h"

#include <phparameter/PHParameters.h>

#include <iostream>  // for operator<<, endl, basic_ostream
#include <limits>

class PHCompositeNode;
class PHG4Detector;

//_______________________________________________________________________
PHG4BeamlineMagnetSubsystem::PHG4BeamlineMagnetSubsystem(const std::string& na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4BeamlineMagnetSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  /* create magnet */
  m_Detector = new PHG4BeamlineMagnetDetector(this, topNode, GetParams(), Name(), GetLayer());

  return 0;
}

//_______________________________________________________________________
int PHG4BeamlineMagnetSubsystem::process_event(PHCompositeNode* /*topNode*/)
{
  return 0;
}

//_______________________________________________________________________
PHG4Detector* PHG4BeamlineMagnetSubsystem::GetDetector() const
{
  return m_Detector;
}

void PHG4BeamlineMagnetSubsystem::SetDefaultParameters()
{
  set_default_string_param("magtype", "");

  set_default_double_param("field_y", std::numeric_limits<double>::quiet_NaN());
  set_default_double_param("fieldgradient", std::numeric_limits<double>::quiet_NaN());

  set_default_double_param("length", 100);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("radius", 100);
  set_default_double_param("thickness", 100);
  set_default_double_param("tmin", std::numeric_limits<double>::quiet_NaN());
  set_default_double_param("tmax", std::numeric_limits<double>::quiet_NaN());

  set_default_int_param("lengthviarapidity", 1);

  set_default_string_param("material", "G4_Galactic");
}

void PHG4BeamlineMagnetSubsystem::Print(const std::string& /*what*/) const
{
  std::cout << Name() << " Parameters: " << std::endl;
  if (!BeginRunExecuted())
  {
    std::cout << "Need to execute BeginRun() before parameter printout is meaningful" << std::endl;
    std::cout << "To do so either run one or more events or on the command line execute: " << std::endl;
    std::cout << "Fun4AllServer *se = Fun4AllServer::instance();" << std::endl;
    std::cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << std::endl;
    std::cout << "g4->InitRun(se->topNode());" << std::endl;
    std::cout << "PHG4BeamlineMagnetSubsystem *cyl = (PHG4BeamlineMagnetSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << std::endl;
    std::cout << "cyl->Print()" << std::endl;
    return;
  }
  GetParams()->Print();
  return;
}
