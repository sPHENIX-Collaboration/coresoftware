#include "PHG4PhenixEventAction.h"
#include "PHG4EventAction.h"

#include <phool/PHTimer.h>    // for PHTimer

#include <iostream>           // for operator<<, endl, basic_ostream, cout

const int VERBOSE = 0;

PHG4PhenixEventAction::PHG4PhenixEventAction() :
  _timer( PHTimeServer::get()->insert_new( "PHG4PhenixEventAction" ) )
{}

PHG4PhenixEventAction::~PHG4PhenixEventAction()
{
  while (actions_.begin() != actions_.end())
    {
      delete actions_.back();
      actions_.pop_back();
    }
}

//_________________________________________________________________
void PHG4PhenixEventAction::BeginOfEventAction( const G4Event* event )
{
  _timer.get()->restart();

  if ( VERBOSE ) std::cout << "PHG4PhenixEventAction::BeginOfEventAction" << std::endl;
  
  // loop over registered actions, and process
  for( ActionList::const_iterator iter = actions_.begin(); iter != actions_.end(); ++iter )
  {
    if(*iter)
    {
      if ( VERBOSE ) std::cout << "PHG4PhenixEventAction::BeginOfEventAction - processing " << *iter << std::endl;
      (*iter)->BeginOfEventAction( event );
    }
  }

}

//_________________________________________________________________
void PHG4PhenixEventAction::EndOfEventAction( const G4Event* event )
{

  if ( VERBOSE ) std::cout << "PHG4PhenixEventAction::EndOfEventAction" << std::endl;

  // loop over registered actions, and process
  for( ActionList::const_iterator iter = actions_.begin(); iter != actions_.end(); ++iter )
  {
    if(*iter)
    {
      if ( VERBOSE ) std::cout << "PHG4PhenixEventAction::EndOfEventAction - processing " << *iter << std::endl;
      (*iter)->EndOfEventAction( event );
    }
  }

  _timer.get()->stop();
}
