#include "PHG4PhenixSteppingAction.h"
#include "PHG4SteppingAction.h"

PHG4PhenixSteppingAction::~PHG4PhenixSteppingAction()
{
  while (actions_.begin() != actions_.end())
    {
      delete actions_.back();
      actions_.pop_back();
    }
}


//_________________________________________________________________
void PHG4PhenixSteppingAction::UserSteppingAction( const G4Step* aStep )
{
  // loop over registered actions, and process
  bool hit_was_used = false;
  for( ActionList::const_iterator iter = actions_.begin(); iter != actions_.end(); ++iter )
  {
    if(*iter)
    {
      hit_was_used |= (*iter)->UserSteppingAction( aStep, hit_was_used );
    }
  }

}
