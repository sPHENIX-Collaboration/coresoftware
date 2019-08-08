// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4MVTX_PHG4VOUTERTRACKERSTEPPINGACTION_H
#define G4MVTX_PHG4VOUTERTRACKERSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4OuterTrackerDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class G4VPhysicalVolume;
class PHParameters;

class PHG4OuterTrackerSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4OuterTrackerSteppingAction(PHG4OuterTrackerDetector *, const int layer);

  //! destructor
  virtual ~PHG4OuterTrackerSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4OuterTrackerDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  PHG4HitContainer *m_AbsorberhitContainer;
  PHG4Hit *m_Hit;
  PHG4Shower *m_SaveShower;
  PHG4HitContainer *m_SaveHitContainer;
  PHG4HitContainer *m_AbsorberHits;

  G4VPhysicalVolume *savevolpre;
  G4VPhysicalVolume *savevolpost;
  int savetrackid;
  int saveprestepstatus;
  int savepoststepstatus;

  int layer_id;
  //const PHParameters *params;
};

#endif  // G4MVTX_PHG4VOUTERTRACKERSTEPPINGACTION_H
