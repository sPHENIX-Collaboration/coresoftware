#include "PHG4ConeSubsystem.h"
#include "PHG4ConeDetector.h"
#include "PHG4ConeSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Subsystem.h>          // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4String.hh>         // for G4String
#include <Geant4/G4SystemOfUnits.hh>  // for cm
#include <Geant4/G4Types.hh>          // for G4double

#include <cmath>  // for tan, atan, exp, M_PI
#include <sstream>

class PHG4Detector;
class PHG4SteppingAction;

using namespace std;

//_______________________________________________________________________
PHG4ConeSubsystem::PHG4ConeSubsystem(const std::string& name, const int lyr)
  : PHG4DetectorSubsystem(name,lyr)
  , detector_(nullptr)
  , layer(lyr)
  , detector_type(name)
{
  InitializeParameters();
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
}

//_______________________________________________________________________
int PHG4ConeSubsystem::InitRunSubsystem(PHCompositeNode* topNode)
{
  // create detector
  detector_ = new PHG4ConeDetector(this, topNode, GetParams(), Name(), layer);
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
    PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));

    string nodename;
    if (SuperDetector() != "NONE")
    {
       PHNodeIterator iter_dst(dstNode);
      PHCompositeNode *superSubNode = dynamic_cast<PHCompositeNode *>(iter_dst.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        dstNode->addNode(superSubNode);
      }
      dstNode = superSubNode;
      PHNodeIterator iter_run(runNode);
      superSubNode = dynamic_cast<PHCompositeNode *>(iter_run.findFirst("PHCompositeNode", SuperDetector()));
      if (!superSubNode)
      {
        superSubNode = new PHCompositeNode(SuperDetector());
        runNode->addNode(superSubNode);
      }
     nodename = "G4HIT_" + SuperDetector();
    }
    else
    {
      nodename = "G4HIT_" + detector_type + "_" + std::to_string(layer);
    }
    // create hit list
    PHG4HitContainer* cone_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename);
    if (!cone_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(cone_hits = new PHG4HitContainer(nodename), nodename, "PHObject"));
    }
    // create stepping action
    m_SteppingAction = new PHG4ConeSteppingAction(detector_);

  }
  return 0;
}

//_______________________________________________________________________
int PHG4ConeSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4ConeSubsystem::GetDetector(void) const
{
  return detector_;
}

//_______________________________________________________________________
void PHG4ConeSubsystem::Set_eta_range(G4double etaMin, G4double etaMax)
{
  G4double thetaMin = 2 * atan(exp(-etaMax));
  G4double thetaMax = 2 * atan(exp(-etaMin));

  G4double z1 = get_double_param("place_z") - get_double_param("length");
  G4double z2 = get_double_param("place_z") + get_double_param("length");

set_double_param("rmin1", z1 * tan(thetaMin));
set_double_param("rmax1", z1 * tan(thetaMax));

set_double_param("rmin2", z2 * tan(thetaMin));
set_double_param("rmax2", z2 * tan(thetaMax));
}

void PHG4ConeSubsystem::SetDefaultParameters()
{
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);

  set_default_double_param("length", NAN);
  set_default_double_param("rmin1", NAN);
  set_default_double_param("rmax1", NAN);
  set_default_double_param("rmin2", NAN);
  set_default_double_param("rmax2", NAN);
  set_default_double_param("sphi", 0.);
  set_default_double_param("dphi", 360.); // degrees
  set_default_double_param("rot_z", 0);

  set_default_string_param("material", "WorldMaterial");
}

void PHG4ConeSubsystem::SetR1(const G4double min, const G4double max)
{
  set_double_param("rmin1", min);
  set_double_param("rmax1", max);
}

void PHG4ConeSubsystem::SetR2(const G4double min, const G4double max)
{
  set_double_param("rmin2", min);
  set_double_param("rmax2", max);
}

void PHG4ConeSubsystem::SetZlength(const G4double a)
{
  set_double_param("length",a);
}

void PHG4ConeSubsystem::SetPhi(const G4double a, const G4double b)
{
  set_double_param("sphi", a);
  set_double_param("dphi", b);
}

void PHG4ConeSubsystem::SetPlaceZ(const G4double dbl)
{
set_double_param("place_z",dbl);
}

  void PHG4ConeSubsystem::SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
set_double_param("place_x",place_x);
set_double_param("place_y",place_y);
set_double_param("place_z",place_z);
  }

void PHG4ConeSubsystem::SetZRot(const G4double dbl)
{
set_double_param("rot_z",dbl);
}

void PHG4ConeSubsystem::SetMaterial(const std::string &mat)
{
set_string_param("material",mat);
}
