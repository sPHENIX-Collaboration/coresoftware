// This is the subsystem code (cc file) for the hcal prototype detector
// created on 1/27/2014, Liang, HeXC
//
#include "PHG4HcalPrototypeSubsystem.h"
#include "PHG4HcalPrototypeDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4HcalPrototypeSteppingAction.h"

#include <g4main/PHG4HitContainer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4HcalPrototypeSubsystem::PHG4HcalPrototypeSubsystem( const std::string &name, const int lyr ):
  PHG4Subsystem( name ),
  detector_( 0 ),
  steppingAction_( nullptr ),
  eventAction_(nullptr),
  place_in_x(0),
  place_in_y(0),
  place_in_z(0),
  rot_in_x(0),
  rot_in_y(0),
  rot_in_z(0),
  material("G4_AIR"),  // default - almost nothing
  active(0),
  absorberactive(0),
  layer(lyr),
  blackhole(0),
  detector_type(name),
  superdetector("NONE")
{

  // put the layer into the name so we get unique names
  // for multiple layers
  ostringstream nam;
  nam << name << "_" << lyr;
  Name(nam.str().c_str());
  for (int i = 0; i < 3; i++)
    {
      dimension[i] = 100.0 * cm;
    }
}

//_______________________________________________________________________
int PHG4HcalPrototypeSubsystem::Init( PHCompositeNode* topNode )
{
  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  detector_ = new PHG4HcalPrototypeDetector(topNode, Name(), layer);
  //  detector_->SetPlace(place_in_x, place_in_y, place_in_z);
  //  detector_->SetXRot(rot_in_x);
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
	  nodename <<  "G4HIT_" << superdetector;
	}
      else
	{
	  nodename <<  "G4HIT_" << detector_type << "_" << layer;
	}
      // create hit list
      PHG4HitContainer* block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
      if ( !block_hits )
	{

	  dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));

	}
      if (absorberactive)
	{
	  nodename.str("");
	  if (superdetector != "NONE")
	    {
	      nodename <<  "G4HIT_ABSORBER_" << superdetector;
	    }
	  else
	    {
	      nodename <<  "G4HIT_ABSORBER_" << detector_type << "_" << layer;
	    }
	}
      block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
      if ( !block_hits )
	{

	  dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));

	}
      // create stepping action
      steppingAction_ = new PHG4HcalPrototypeSteppingAction(detector_);

      eventAction_ = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
    }
  if (blackhole && !active)
    {
      steppingAction_ = new PHG4HcalPrototypeSteppingAction(detector_);
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4HcalPrototypeSubsystem::process_event( PHCompositeNode * topNode )
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
    if (steppingAction_)
        {
            steppingAction_->SetInterfacePointers( topNode );
        }
    return 0;
}


//_______________________________________________________________________
PHG4Detector* PHG4HcalPrototypeSubsystem::GetDetector( void ) const
{
    return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4HcalPrototypeSubsystem::GetSteppingAction( void ) const
{
    return steppingAction_;
}

