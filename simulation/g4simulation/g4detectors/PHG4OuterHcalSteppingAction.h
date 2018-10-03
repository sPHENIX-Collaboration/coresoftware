// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4DETECTORS_PHG4OUTERHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4OUTERHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4VPhysicalVolume;
class PHG4OuterHcalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4OuterHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4OuterHcalSteppingAction(PHG4OuterHcalDetector *, const PHParameters *parameters);

  //! destructor
  virtual ~PHG4OuterHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  virtual int Init();

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

  double GetLightCorrection(const double r) const;

  void FieldChecker(const G4Step *);
  void EnableFieldChecker(const int i = 1) { m_EnableFieldCheckerFlag = i; }
 
private:
  //! pointer to the detector
  PHG4OuterHcalDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *m_Hits;
  PHG4HitContainer *m_AbsorberHits;
  PHG4Hit *m_Hit;
  const PHParameters *m_Params;
  PHG4HitContainer *m_SaveHitContainer;
  PHG4Shower *m_SaveShower;
  G4VPhysicalVolume *m_SaveVolPre;
  G4VPhysicalVolume *m_SaveVolPost;
  int m_SaveTrackId;
  int m_SavePreStepStatus;
  int m_SavePostStepStatus;
  int m_EnableFieldCheckerFlag;

  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActiveFlag;
  int m_IsBlackHoleFlag;
  int m_NScintiPlates;
  int m_LightScintModelFlag;

  double m_LightBalanceInnerCorr;
  double m_LightBalanceInnerRadius;
  double m_LightBalanceOuterCorr;
  double m_LightBalanceOuterRadius;
};

#endif  // PHG4OuterHcalSteppingAction_h
