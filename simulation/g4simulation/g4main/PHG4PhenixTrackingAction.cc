#include "PHG4PhenixTrackingAction.h"
#include "PHG4TrackingAction.h"

#include <iostream>

void
PHG4PhenixTrackingAction::PreUserTrackingAction( const G4Track* track )
{

  if ( Verbosity()>0 ) std::cout << "PHG4PhenixTrackingAction::PreUserTrackingAction" << std::endl;
  
  // loop over registered actions, and process
  for( ActionList::const_iterator iter = actions_.begin(); iter != actions_.end(); iter++ )
    {
      if (*iter)
	{
	  if ( Verbosity()>0 ) std::cout << "PHG4PhenixTrackingAction::PreUserTrackingAction - processing " << *iter << std::endl;
	  (*iter)->PreUserTrackingAction(track);
	}
    }

}

void
PHG4PhenixTrackingAction::PostUserTrackingAction( const G4Track* track )
{

  if ( Verbosity()>0 ) std::cout << "PHG4PhenixTrackingAction::PostUserTrackingAction" << std::endl;

  // loop over registered actions, and process
  for( ActionList::const_iterator iter = actions_.begin(); iter != actions_.end(); iter++ )
    {
      if (*iter)
	{
	  if ( Verbosity()>0 ) std::cout << "PHG4PhenixTrackingAction::PostUserTrackingAction - processing " << *iter << std::endl;
	  (*iter)->PostUserTrackingAction(track);
	}
    }

}
