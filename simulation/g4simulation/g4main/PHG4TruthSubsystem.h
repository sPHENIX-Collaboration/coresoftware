#ifndef PHG4TruthSubsystem_h
#define PHG4TruthSubsystem_h

#include "PHG4Subsystem.h"
#include <string>

class PHG4TruthSteppingAction;
class PHG4TruthTrackingAction;
class PHG4TruthEventAction;

class PHG4TruthSubsystem: public PHG4Subsystem
{

  public:

  //! constructor
  PHG4TruthSubsystem( const std::string &name = "TRUTH" );

  //! destructor
  virtual ~PHG4TruthSubsystem( void )
  {}

  //! init
  int InitRun(PHCompositeNode *);

  //! event processing
  int process_event(PHCompositeNode *);

  //! Clean up after each event.
  int ResetEvent(PHCompositeNode *);

  //! accessors (reimplemented)
  virtual PHG4EventAction* GetEventAction( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const;
  virtual PHG4TrackingAction* GetTrackingAction( void ) const;

  private:

  PHG4TruthEventAction* eventAction_;
  PHG4TruthSteppingAction* steppingAction_;
  PHG4TruthTrackingAction* trackingAction_;

};

#endif
