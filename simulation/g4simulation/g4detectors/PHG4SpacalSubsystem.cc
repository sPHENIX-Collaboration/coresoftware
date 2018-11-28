/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */
#include "PHG4SpacalSubsystem.h"

#include "PHG4CylinderGeom.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4FullProjSpacalDetector.h"
#include "PHG4FullProjTiltedSpacalDetector.h"
//#include "PHG4ProjSpacalDetector.h"
#include "PHG4SpacalDetector.h"
#include "PHG4SpacalSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4SpacalSubsystem::PHG4SpacalSubsystem(const std::string& na, const int lyr)
  : PHG4DetectorSubsystem(na, lyr)
  , detector_(nullptr)
  , steppingAction_(nullptr)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4SpacalSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  // create hit list only for active layers
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  // create detector
  //  _geom.set_layer( layer );
  //  if (lengthViaRapidityCoverage)
  //    {
  //       const double half_length =  PHG4Utils::GetLengthForRapidityCoverage(_geom.get_max_radius());
  //       _geom.set_zmin(-half_length);
  //       _geom.set_zmax(+half_length);
  //    }

  switch (GetParams()->get_int_param("config"))
  {
  case PHG4CylinderGeom_Spacalv1::kNonProjective:
    if (Verbosity() > 0) cout << "PHG4SpacalSubsystem::InitRun - use PHG4SpacalDetector" << endl;
    detector_ = new PHG4SpacalDetector(topNode, Name(), GetParams(), GetLayer());
    break;

  case PHG4CylinderGeom_Spacalv1::kProjective_PolarTaper:
    cout << "PHG4SpacalSubsystem::InitRun - PHG4ProjSpacalDetector is obsolete" << endl;
    exit(10);
//    detector_ = new PHG4ProjSpacalDetector(topNode, Name(), GetParams(), GetLayer());
    break;

  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper:
  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper_SameLengthFiberPerTower:
    if (Verbosity() > 0) cout << "PHG4SpacalSubsystem::InitRun - use PHG4FullProjSpacalDetector" << endl;
    detector_ = new PHG4FullProjSpacalDetector(topNode, Name(), GetParams(), GetLayer());
    break;

  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper_Tilted:
  case PHG4CylinderGeom_Spacalv1::kFullProjective_2DTaper_Tilted_SameLengthFiberPerTower:
    if (Verbosity() > 0) cout << "PHG4SpacalSubsystem::InitRun - use PHG4FullProjTiltedSpacalDetector" << endl;
    detector_ = new PHG4FullProjTiltedSpacalDetector(topNode, Name(), GetParams(), GetLayer());
    break;

  default:
    cout << "PHG4SpacalSubsystem::InitRun - unknown option exiting" << endl;
    exit(1);
    break;
  }

  detector_->SetActive(GetParams()->get_int_param("active"));
  detector_->SetAbsorberActive(GetParams()->get_int_param("absorberactive"));
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());

  if (GetParams()->get_int_param("active"))
  {
    ostringstream nodename;
    if (SuperDetector()  != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector() ;
    }
    else
    {
      nodename << "G4HIT_" <<  Name() << "_" << GetLayer();
    }
    PHG4HitContainer* cylinder_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!cylinder_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(cylinder_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    cylinder_hits->AddLayer(GetLayer());
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodename.str("");
      if (SuperDetector() != "NONE")
      {
        nodename << "G4HIT_ABSORBER_" << SuperDetector();
      }
      else
      {
        nodename << "G4HIT_ABSORBER_" << Name() << "_" << GetLayer();
      }
      PHG4HitContainer* cylinder_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
      if (!cylinder_hits)
      {
        dstNode->addNode(new PHIODataNode<PHObject>(cylinder_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
      }
      cylinder_hits->AddLayer(GetLayer());
    }
    steppingAction_ = new PHG4SpacalSteppingAction(detector_);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4SpacalSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4SpacalSubsystem::GetDetector(void) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4SpacalSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
}

void PHG4SpacalSubsystem::Print(const std::string& what) const
{
  detector_->Print(what);
  return;
}

void PHG4SpacalSubsystem::SetDefaultParameters()
{
  set_default_double_param("xpos", 0.);  // translation in 3D
  set_default_double_param("ypos", 0.);  // translation in 3D
  set_default_double_param("zpos", 0.);  // translation in 3D

  set_default_double_param("thickness", 21.00000);
  set_default_double_param("radius", 90.);
  set_default_double_param("zmin", -149.470000);
  set_default_double_param("zmax", 149.470000);
  set_default_int_param("azimuthal_n_sec", 256);

  set_default_int_param("construction_verbose", 0.);
  set_default_int_param("azimuthal_seg_visible", 0.);
  set_default_int_param("virualize_fiber", 0.);
  set_default_int_param("config", static_cast<int>(PHG4CylinderGeom_Spacalv1::kNonProjective));

  set_default_double_param("divider_width", 0); // radial size of the divider between blocks. <=0 means no dividers
  set_default_string_param("divider_mat", "G4_AIR"); // materials of the divider. G4_AIR is equivalent to not installing one in the term of material distribution

  return;
}
