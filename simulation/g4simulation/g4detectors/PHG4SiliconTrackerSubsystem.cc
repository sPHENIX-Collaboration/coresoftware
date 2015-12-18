#include "PHG4SiliconTrackerSubsystem.h"
#include "PHG4SiliconTrackerDetector.h"
#include "PHG4EventActionClearZeroEdep.h"
#include "PHG4SiliconTrackerSteppingAction.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4SiliconTrackerSubsystem::PHG4SiliconTrackerSubsystem( const std::string &name, const int lyr ):
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
  material("G4_AIR"),  // default - almost nothing
  active(0),
  absorberactive(0),
  layer(lyr),
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
int PHG4SiliconTrackerSubsystem::Init( PHCompositeNode* topNode )
{
  if(verbosity>0)
    cout << "PHG4SiliconTrackerSubsystem::Init started" << endl;

  PHNodeIterator iter( topNode );
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST" ));

  // create detector
  cout << "    create detector for layer " << layer << endl;
  cout << "strip tilt " << strip_tilt << " rad " << rad << " strip_tilt * rad " << strip_tilt * rad << endl; 
  detector_ = new PHG4SiliconTrackerDetector(topNode, Name(), layer);
  detector_->Verbosity(verbosity);
  detector_->set_nominal_layer_radius(layer_nominal_radius);
  detector_->set_radius_stagger(radius_stagger);
  detector_->set_N_staggers(N_staggers);
  detector_->set_N_strips_in_sensor_phi(N_strips_per_column);
  detector_->set_strip_tilt(strip_tilt);
  detector_->set_option_double_layer(option_double_layer);
  detector_->set_add_lower_roc(add_lower_roc);
  detector_->SetPlace(place_in_x, place_in_y, place_in_z);
  detector_->SetXRot(rot_in_x);
  detector_->SetYRot(rot_in_y);
  detector_->SetZRot(rot_in_z);
  detector_->SetActive(active);
  detector_->SetAbsorberActive(absorberactive);
  detector_->BlackHole(blackhole);
  detector_->SuperDetector(superdetector);
  detector_->Detector(detector_type);
  detector_->OverlapCheck(overlapcheck);
  cout << "    ------ created detector for " << layer  << endl;

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
      steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_);
    }
  if (blackhole && !active)
    {
      steppingAction_ = new PHG4SiliconTrackerSteppingAction(detector_);
    }
  return 0;

}

//_______________________________________________________________________
int
PHG4SiliconTrackerSubsystem::process_event( PHCompositeNode * topNode )
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
PHG4Detector* PHG4SiliconTrackerSubsystem::GetDetector( void ) const
{
    return detector_;
}

//_______________________________________________________________________
PHG4SteppingAction* PHG4SiliconTrackerSubsystem::GetSteppingAction( void ) const
{
    return steppingAction_;
}

