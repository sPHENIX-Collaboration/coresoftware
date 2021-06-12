#include "PHG4PhenixStackingAction.h"
#include "PHG4StackingAction.h"

PHG4PhenixStackingAction::~PHG4PhenixStackingAction()
{
  while (actions_.begin() != actions_.end())
    {
      delete actions_.back();
      actions_.pop_back();
    }
}


//_________________________________________________________________
G4ClassificationOfNewTrack PHG4PhenixStackingAction::ClassifyNewTrack( const G4Track* aTrack )
{
  // loop over registered actions, and process
  G4ClassificationOfNewTrack retcode = fUrgent;
  for( ActionList::const_iterator iter = actions_.begin(); iter != actions_.end(); ++iter )
  {
    if(*iter)
    {
      G4ClassificationOfNewTrack cret = (*iter)->ClassifyNewTrack(aTrack);
      if (cret != fUrgent)
      {
	retcode = cret;
      }
    }
  }
  return retcode;
}

//_________________________________________________________________
void PHG4PhenixStackingAction::PrepareNewEvent()
{
  // loop over registered actions, and process
  for( ActionList::const_iterator iter = actions_.begin(); iter != actions_.end(); ++iter )
  {
    if(*iter)
    {
      (*iter)->PrepareNewEvent();
    }
  }
  return;
}
