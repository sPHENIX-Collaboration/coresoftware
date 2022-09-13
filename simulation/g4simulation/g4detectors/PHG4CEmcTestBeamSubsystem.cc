#include "PHG4CEmcTestBeamSubsystem.h"
#include "PHG4CEmcTestBeamDetector.h"

#include "PHG4CEmcTestBeamSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Subsystem.h>  // for PHG4Subsystem

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>  // for G4double

#include <sstream>

class PHG4Detector;
class PHG4SteppingAction;

using namespace std;

//_______________________________________________________________________
PHG4CEmcTestBeamSubsystem::PHG4CEmcTestBeamSubsystem(const std::string& name, const int lyr)
  : PHG4Subsystem(name)
  , detector_(nullptr)
  , steppingAction_(nullptr)
  , eventAction_(nullptr)
  , place_in_x(0)
  , place_in_y(0)
  , place_in_z(0)
  , rot_in_x(0)
  , rot_in_y(0)
  , rot_in_z(0)
  , active(0)
  , absorberactive(0)
  , layer(lyr)
  , blackhole(0)
  , detector_type(name)
  , superdetector("NONE")
{
  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str());
  for (double& i : dimension)
  {
    i = 100.0 * cm;
  }
}

//_______________________________________________________________________
int PHG4CEmcTestBeamSubsystem::Init(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4CEmcTestBeamDetector(this, topNode, Name(), layer);
  detector_->SetPlace(place_in_x, place_in_y, place_in_z);
  detector_->SetXRot(rot_in_x);
  detector_->SetYRot(rot_in_y);
  detector_->SetZRot(rot_in_z);
  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->BlackHole(blackhole);
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
    PHG4HitContainer* block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!block_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    if (absorberactive)
    {
      nodename.str("");
      if (superdetector != "NONE")
      {
        nodename << "G4HIT_ABSORBER_" << superdetector;
      }
      else
      {
        nodename << "G4HIT_ABSORBER_" << detector_type << "_" << layer;
      }
    }
    block_hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str());
    if (!block_hits)
    {
      dstNode->addNode(new PHIODataNode<PHObject>(block_hits = new PHG4HitContainer(nodename.str()), nodename.str(), "PHObject"));
    }
    // create stepping action
    steppingAction_ = new PHG4CEmcTestBeamSteppingAction(detector_);

    eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
  }
  if (blackhole && !active)
  {
    steppingAction_ = new PHG4CEmcTestBeamSteppingAction(detector_);
  }
  return 0;
}

//_______________________________________________________________________
int PHG4CEmcTestBeamSubsystem::process_event(PHCompositeNode* topNode)
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
PHG4Detector* PHG4CEmcTestBeamSubsystem::GetDetector() const
{
  return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4CEmcTestBeamSubsystem::GetSteppingAction() const
{
  return steppingAction_;
}

void PHG4CEmcTestBeamSubsystem::SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
{
  place_in_x = place_x * cm;
  place_in_y = place_y * cm;
  place_in_z = place_z * cm;
}

void PHG4CEmcTestBeamSubsystem::SetPlaceZ(const G4double dbl)
{
  place_in_z = dbl * cm;
}

void PHG4CEmcTestBeamSubsystem::SetXRot(const G4double dbl)
{
  rot_in_x = dbl * deg;
}

void PHG4CEmcTestBeamSubsystem::SetYRot(const G4double dbl)
{
  rot_in_y = dbl * deg;
}

void PHG4CEmcTestBeamSubsystem::SetZRot(const G4double dbl)
{
  rot_in_z = dbl * deg;
}
