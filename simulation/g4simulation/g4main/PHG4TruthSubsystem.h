#ifndef G4MAIN_PHG4TRUTHSUBSYSTEM_H
#define G4MAIN_PHG4TRUTHSUBSYSTEM_H

#include "PHG4Subsystem.h"

#include <string>

class PHG4TruthTrackingAction;
class PHG4TruthEventAction;

class PHG4TruthSubsystem : public PHG4Subsystem
{
 public:
  //! constructor
  PHG4TruthSubsystem(const std::string &name = "TRUTH");

  //! destructor
  virtual ~PHG4TruthSubsystem(void)
  {
  }

  //! init
  int InitRun(PHCompositeNode *);

  //! event processing
  int process_event(PHCompositeNode *);

  //! event processing
  virtual int process_after_geant(PHCompositeNode *);

  //! Clean up after each event.
  int ResetEvent(PHCompositeNode *);

  //! accessors (reimplemented)
  virtual PHG4EventAction *GetEventAction(void) const;
  virtual PHG4TrackingAction *GetTrackingAction(void) const;

  //! only save the G4 truth information that is associated with the embedded particle
  void SetSaveOnlyEmbeded(bool b = true) { m_SaveOnlyEmbededFlag = b; };

 private:
  PHG4TruthEventAction *m_EventAction;

  PHG4TruthTrackingAction *m_TrackingAction;

  //! only save the G4 truth information that is associated with the embedded particle
  bool m_SaveOnlyEmbededFlag;
};

#endif
