#include "PHG4SpacalPrototypeSubsystem.h"


#include "PHG4SpacalPrototypeDetector.h"
#include "PHG4FullProjSpacalDetector.h"
#include "PHG4CylinderGeom.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4SpacalPrototypeSteppingAction.h"

#include <phparameter/PHParameters.h>
#include <phparameter/PHParametersContainer.h>

#include <g4main/PHG4Utils.h>
#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <Geant4/globals.hh>

#include <sstream>
#include <cassert>

using namespace std;

//_______________________________________________________________________
PHG4SpacalPrototypeSubsystem::PHG4SpacalPrototypeSubsystem(const std::string &na) :
    PHG4DetectorSubsystem(na,0),
    detector_(nullptr), 
    steppingAction_(nullptr)
{
  InitializeParameters();
}

//_______________________________________________________________________
int
PHG4SpacalPrototypeSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  if (Verbosity() > 0)
    cout
        << "PHG4SpacalPrototypeSubsystem::InitRun - use PHG4SpacalPrototypeDetector"
        << endl;
  detector_ = new PHG4SpacalPrototypeDetector(topNode, GetParams(), Name());

  detector_->SetActive(GetParams()->get_int_param("active"));
  detector_->SetAbsorberActive(GetParams()->get_int_param("absorberactive"));
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
    {
      ostringstream nodename;
      if (SuperDetector() != "NONE")
        {
          nodename << "G4HIT_" << SuperDetector();
        }
      else
        {
          nodename << "G4HIT_" << Name();
        }
      PHG4HitContainer* cylinder_hits = findNode::getClass<PHG4HitContainer>(
          topNode, nodename.str().c_str());
      if (!cylinder_hits)
        {
          dstNode->addNode(
              new PHIODataNode<PHObject>(
                  cylinder_hits = new PHG4HitContainer(nodename.str()),
                  nodename.str().c_str(), "PHObject"));
        }
      cylinder_hits->AddLayer(0);
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
          PHG4HitContainer* cylinder_hits =
              findNode::getClass<PHG4HitContainer>(topNode,
                  nodename.str().c_str());
          if (!cylinder_hits)
            {
              dstNode->addNode(
                  new PHIODataNode<PHObject>(cylinder_hits =
                      new PHG4HitContainer(nodename.str()),
                      nodename.str().c_str(), "PHObject"));
            }
          cylinder_hits->AddLayer(0);
        }
      steppingAction_ = new PHG4SpacalPrototypeSteppingAction(detector_);
    }

  return 0;

}

//_______________________________________________________________________
int
PHG4SpacalPrototypeSubsystem::process_event(PHCompositeNode* topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
    {
      steppingAction_->SetInterfacePointers(topNode);
    }
  return 0;

}

//_______________________________________________________________________
PHG4Detector*
PHG4SpacalPrototypeSubsystem::GetDetector(void) const
{
  return detector_;
}

void
PHG4SpacalPrototypeSubsystem::Print(const std::string &what) const
{
  detector_->Print(what);
  cout << Name() << " Parameters: " << endl;
  if (! BeginRunExecuted())
    {
      cout << "Need to execute BeginRun() before parameter printout is meaningful" << endl;
      cout << "To do so either run one or more events or on the command line execute: " << endl;
      cout << "Fun4AllServer *se = Fun4AllServer::instance();" << endl;
      cout << "PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");" << endl;
      cout << "g4->InitRun(se->topNode());" << endl;
      cout << "PHG4SpacalPrototypeSubsystem *sys = (PHG4SpacalPrototypeSubsystem *) g4->getSubsystem(\"" << Name() << "\");" << endl;
      cout << "sys->Print()" << endl;
      return;
    }
  GetParams()->Print();
  return;
}

void
PHG4SpacalPrototypeSubsystem::SetDefaultParameters()
{
  set_default_double_param("xpos", 0.); // translation in 3D
  set_default_double_param("ypos", 0.); // translation in 3D
  set_default_double_param("zpos", 0.); // translation in 3D
  set_default_double_param("z_rotation_degree", 0.); // roation in the vertical plane
  set_default_int_param("construction_verbose", 0.); // roation in the vertical plane
  return;
}
