#include "PHG4ConeSubsystem.h"
#include "PHG4ConeDetector.h"
#include "PHG4ConeSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"

#include <g4main/PHG4HitContainer.h>

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
  : PHG4Subsystem(name)
  , detector_(0)
  , steppingAction_(nullptr)
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
  , active(0)
  , layer(lyr)
  , detector_type(name)
  , superdetector("NONE")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
}

//_______________________________________________________________________
int PHG4ConeSubsystem::Init(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4ConeDetector(this, topNode, Name(), layer);
  detector_->SetR1(rMin1, rMax1);
  detector_->SetR2(rMin2, rMax2);
  detector_->SetZlength(dZ);
  detector_->SetPhi(sPhi, dPhi);
  detector_->SetPlace(place_in_x, place_in_y, place_in_z);
  detector_->SetZRot(rot_in_z);
  detector_->SetMaterial(material);
  detector_->SetActive(active);
  detector_->SuperDetector(superdetector);
  detector_->OverlapCheck(CheckOverlap());
  if (active)
  {
    ostringstream nodename;
    if (superdetector != "NONE")
    {
      nodename << "G4HIT_" << superdetector;
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
    steppingAction_ = new PHG4ConeSteppingAction(detector_);

    eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
  }
  return 0;
}

//_______________________________________________________________________
int PHG4ConeSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4ConeSubsystem::GetDetector(void) const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4ConeSubsystem::GetSteppingAction(void) const
{
  return steppingAction_;
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

  rMin2 = z2 * tan(thetaMin);
  rMax2 = z2 * tan(thetaMax);
}
