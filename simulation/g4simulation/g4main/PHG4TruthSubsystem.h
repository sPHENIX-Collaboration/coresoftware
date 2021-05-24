// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4TRUTHSUBSYSTEM_H
#define G4MAIN_PHG4TRUTHSUBSYSTEM_H

#include "PHG4Subsystem.h"

#include <string>

class PHCompositeNode;
class PHG4EventAction;
class PHG4TrackingAction;
class PHG4TruthTrackingAction;
class PHG4TruthEventAction;

class PHG4TruthSubsystem : public PHG4Subsystem
{
 public:
  //! constructor
  PHG4TruthSubsystem(const std::string &name = "TRUTH");

  //! destructor
  ~PHG4TruthSubsystem(void) override
  {
  }

  //! init
  int InitRun(PHCompositeNode *) override;

  //! event processing
  int process_event(PHCompositeNode *) override;

  //! event processing
  int process_after_geant(PHCompositeNode *) override;

  //! Clean up after each event.
  int ResetEvent(PHCompositeNode *) override;

  //! accessors (reimplemented)
  PHG4EventAction *GetEventAction(void) const override;
  PHG4TrackingAction *GetTrackingAction(void) const override;

  //! only save the G4 truth information that is associated with the embedded particle
  void SetSaveOnlyEmbeded(bool b = true) { m_SaveOnlyEmbededFlag = b; };

 private:
  PHG4TruthEventAction *m_EventAction;

  PHG4TruthTrackingAction *m_TrackingAction;

  //! only save the G4 truth information that is associated with the embedded particle
  bool m_SaveOnlyEmbededFlag;
};

#endif
