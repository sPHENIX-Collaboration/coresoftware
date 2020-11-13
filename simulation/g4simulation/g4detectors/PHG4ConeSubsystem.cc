#include "PHG4ConeSubsystem.h"
#include "PHG4ConeDetector.h"
#include "PHG4ConeSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"

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
  , eventAction_(nullptr)
  , place_in_x(0)
  , place_in_y(0)
  , place_in_z(0)
  , rot_in_z(0)
  , rMin1(5 * cm)
  , rMax1(100 * cm)
  , rMin2(5 * cm)
  , rMax2(200 * cm)
  , dZ(100 * cm)
  , sPhi(0)
  , dPhi(2 * M_PI)
  , material("Silicon")
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
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4ConeDetector(this, topNode, GetParams(), Name(), layer);
  detector_->SetR1(rMin1, rMax1);
  detector_->SetR2(rMin2, rMax2);
  detector_->SetZlength(dZ);
  detector_->SetPhi(sPhi, dPhi);
  detector_->SetPlace(place_in_x, place_in_y, place_in_z);
  detector_->SetZRot(rot_in_z);
  detector_->SetMaterial(material);
  detector_->SetActive(GetParams()->get_int_param("active"));
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
      nodename << "G4HIT_" << detector_type << "_" << layer;
    }
    // create hit list
    PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
    if (!block_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject"));
    }
    // create stepping action
    m_SteppingAction = new PHG4ConeSteppingAction(detector_);

    eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
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

  G4double z1 = place_in_z - dZ;
  G4double z2 = place_in_z + dZ;

  rMin1 = z1 * tan(thetaMin);
  rMax1 = z1 * tan(thetaMax);
set_double_param("rmin1", rMin1);
set_double_param("rmax1", rMax1);

  rMin2 = z2 * tan(thetaMin);
  rMax2 = z2 * tan(thetaMax);
set_double_param("rmin2", rMin2);
set_double_param("rmax2", rMax2);
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
  set_default_double_param("sphi", NAN);
  set_default_double_param("dphi", NAN);
  set_default_double_param("rot_z", NAN);

  set_default_string_param("material", "WorldMaterial");
}

void PHG4ConeSubsystem::SetR1(const G4double min, const G4double max)
{
  rMin1 = min; 
  rMax1=max;
  set_double_param("rmin1", min);
  set_double_param("rmax1", max);
}

void PHG4ConeSubsystem::SetR2(const G4double min, const G4double max)
{
  rMin2 = min; 
  rMax2=max;
  set_double_param("rmin2", min);
  set_double_param("rmax2", max);
}

void PHG4ConeSubsystem::SetZlength(const G4double a)
{
  dZ = a;
  set_double_param("length",a);
}

void PHG4ConeSubsystem::SetPhi(const G4double a, const G4double b)
{
sPhi = a; 
dPhi = b;
  set_double_param("sphi", a);
  set_double_param("dphi", b);
}

void PHG4ConeSubsystem::SetPlaceZ(const G4double dbl)
{
place_in_z = dbl;
set_double_param("place_z",dbl);
}

  void PHG4ConeSubsystem::SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
    place_in_x = place_x;
    place_in_y = place_y;
    place_in_z = place_z;
set_double_param("place_x",place_x);
set_double_param("place_y",place_y);
set_double_param("place_z",place_z);
  }

void PHG4ConeSubsystem::SetZRot(const G4double dbl)
{
rot_in_z = dbl;
set_double_param("rot_z",dbl);
}

void PHG4ConeSubsystem::SetMaterial(const std::string &mat)
{
material = mat;
set_string_param("material",mat);
}
