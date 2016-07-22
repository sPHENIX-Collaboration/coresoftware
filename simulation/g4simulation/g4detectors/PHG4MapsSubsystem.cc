#include "PHG4MapsSubsystem.h"
#include "PHG4MapsDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4MapsSteppingAction.h"

#include "Geant4/G4GDMLParser.hh"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4MapsSubsystem::PHG4MapsSubsystem( const std::string &name, const int lyr, int in_stave_type ):
  PHG4Subsystem( name ),
  detector_( 0 ),
  steppingAction_( NULL ),
  eventAction_(NULL),
  place_in_x(0),
  place_in_y(0),
  place_in_z(0),
  rot_in_x(0),
  rot_in_y(0),
  rot_in_z(0),
  layer(lyr),
  stave_type(in_stave_type),
  material("G4_AIR"),  // default - almost nothing
  active(0),
  absorberactive(0),
  blackhole(0),
  detector_type(name),
  superdetector(name)
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
int PHG4MapsSubsystem::Init( PHCompositeNode* topNode )
{
  if(verbosity>0)
    cout << "PHG4MapsSubsystem::Init started" << endl;

  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  // These values are set from the calling macro using the setters defined in the .h file
  cout << "    create MAPS detector for layer " << layer << endl;
  detector_ = new PHG4MapsDetector(topNode, Name(), layer, stave_type);
  detector_->Verbosity(2);
  detector_->set_nominal_layer_radius(layer_nominal_radius);
  detector_->set_pixel_x(pixel_x);
  detector_->set_pixel_y(pixel_y);
  //detector_->set_stave_type(stave_type);
  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->BlackHole(blackhole);
  detector_->SuperDetector(superdetector);
  detector_->Detector(detector_type);
  detector_->OverlapCheck(overlapcheck);
  cout << "    ------ created detector for " << layer  << " name " << Name() << endl;

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
      cout << PHWHERE << "creating hits node " << nodename.str() << endl;

      PHG4EventActionClearZeroEdep *eventaction = new PHG4EventActionClearZeroEdep(topNode, nodename.str());
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
	  block_hits =  findNode::getClass<PHG4HitContainer>( topNode , nodename.str().c_str());
	  if ( !block_hits )
	    {

	      dstNode->addNode( new PHIODataNode<PHObject>( block_hits = new PHG4HitContainer(nodename.str()), nodename.str().c_str(), "PHObject" ));

	    }
	  eventaction->AddNode(nodename.str());
	}
      eventAction_ = dynamic_cast<PHG4EventAction *> (eventaction);
      // create stepping action
      steppingAction_ = new PHG4MapsSteppingAction(detector_);
    }
  if (blackhole && !active)
    {
      steppingAction_ = new PHG4MapsSteppingAction(detector_);
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4MapsSubsystem::process_event( PHCompositeNode * topNode )
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
PHG4Detector* PHG4MapsSubsystem::GetDetector( void ) const
{
    return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4MapsSubsystem::GetSteppingAction( void ) const
{
    return steppingAction_;
}

