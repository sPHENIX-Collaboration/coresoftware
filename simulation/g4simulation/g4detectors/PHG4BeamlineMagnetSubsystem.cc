#include "PHG4BeamlineMagnetSubsystem.h"
#include "PHG4BeamlineMagnetDetector.h"

#include <phparameter/PHParameters.h>

#include <cmath>     // for NAN
#include <iostream>  // for operator<<, endl, basic_ostream

class PHCompositeNode;
class PHG4Detector;

using namespace std;

//_______________________________________________________________________
PHG4BeamlineMagnetSubsystem::PHG4BeamlineMagnetSubsystem(const std::string& na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
  , detector_(nullptr)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4BeamlineMagnetSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  /* create magnet */
  detector_ = new PHG4BeamlineMagnetDetector(this, topNode, GetParams(), Name(), GetLayer());

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
  return detector_;
}

void PHG4BeamlineMagnetSubsystem::SetDefaultParameters()
{
  set_default_string_param("magtype", "");

  set_default_double_param("field_y", NAN);
  set_default_double_param("fieldgradient", NAN);

  set_default_double_param("length", 100);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("radius", 100);
  set_default_double_param("thickness", 100);
  set_default_double_param("tmin", NAN);
  set_default_double_param("tmax", NAN);

  set_default_int_param("lengthviarapidity", 1);

  set_default_string_param("material", "G4_Galactic");
}

void PHG4BeamlineMagnetSubsystem::Print(const string& /*what*/) const
{
  cout << Name() << " Parameters: " << endl;
  if (!BeginRunExecuted())
  {
    cout << "Need to execute BeginRun() before parameter printout is meaningful" << endl;
    cout << "To do so either run one or more events or on the command line execute: " << endl;
    cout << "Fun4AllServer *se = Fun4AllServer::instance();" << endl;
    cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << endl;
    cout << "g4->InitRun(se->topNode());" << endl;
    cout << "PHG4BeamlineMagnetSubsystem *cyl = (PHG4BeamlineMagnetSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << endl;
    cout << "cyl->Print()" << endl;
    return;
  }
  GetParams()->Print();
  return;
}
