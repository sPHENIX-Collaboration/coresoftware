#include "PHG4FPbScSubsystem.h"
#include "PHG4FPbScDetector.h"
#include "PHG4FPbScSteppingAction.h"
#include "PHG4FPbScRegionSteppingAction.h"
#include "PHG4EventActionClearZeroEdep.h"

#include <g4main/PHG4HitContainer.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4FPbScSubsystem::PHG4FPbScSubsystem( const string &name ):
PHG4Subsystem( name ),
detector_( NULL ), steppingAction_(NULL), eventAction_(NULL), x_position(0), y_position(0), z_position(0)
{
}

//_______________________________________________________________________
int PHG4FPbScSubsystem::Init( PHCompositeNode* topNode )
{
 
  std::ostringstream hitnodename;
  hitnodename <<  "G4HIT_" << Name(); 
 
  // create hit list
  PHG4HitContainer* fcal_hits =  findNode::getClass<PHG4HitContainer>( topNode , ("G4HIT_"+Name()).c_str() );
  if ( !fcal_hits )
  {
    
    PHNodeIterator iter( topNode );
    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST" ));
    dstNode->addNode( new PHIODataNode<PHObject>( fcal_hits = new PHG4HitContainer(("G4HIT_"+Name())), ("G4HIT_"+Name()),"PHObject" ));
  }
  
  // create detector
  detector_ = new PHG4FPbScDetector(topNode, Name());
  detector_->set_Place(x_position, y_position, z_position);
  detector_->Verbosity(Verbosity());
  
  // create stepping action
  steppingAction_ = new PHG4FPbScSteppingAction(detector_);

  eventAction_ = new PHG4EventActionClearZeroEdep(topNode, hitnodename.str());

  return 0;
  
}

//_______________________________________________________________________
int PHG4FPbScSubsystem::process_event( PHCompositeNode* topNode )
{
  if ( PHG4FPbScRegionSteppingAction* p = dynamic_cast<PHG4FPbScRegionSteppingAction*>(detector_->GetSteppingAction()) )
    p->SetInterfacePointers( topNode );
  else
    steppingAction_->SetInterfacePointers( topNode );

  return 0;
  
}

//_______________________________________________________________________
PHG4Detector* PHG4FPbScSubsystem::GetDetector( void ) const
{ return detector_; }

//_______________________________________________________________________
PHG4SteppingAction* PHG4FPbScSubsystem::GetSteppingAction( void ) const
{ return steppingAction_; }

