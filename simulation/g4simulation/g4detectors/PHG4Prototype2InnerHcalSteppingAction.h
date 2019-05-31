// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROTOTYPE2INNERHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4PROTOTYPE2INNERHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class PHG4Hit;
class PHG4HitContainer;
class PHParameters;
class PHG4Prototype2InnerHcalDetector;
class PHG4Shower;

class PHG4Prototype2InnerHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4Prototype2InnerHcalSteppingAction(PHG4Prototype2InnerHcalDetector *, const PHParameters *parameters);

  //! dtor
  virtual ~PHG4Prototype2InnerHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

  double GetLightCorrection(const double r) const;

 private:
  //! pointer to the detector
  PHG4Prototype2InnerHcalDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer;
  PHG4HitContainer *m_AbsorberHitContainer;
  PHG4Hit *m_Hit;
  const PHParameters *m_Params;
  PHG4HitContainer *m_SaveHitContainer;
  PHG4Shower *m_SaveShower;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_AbsorberTruthFlag;
  int m_IsActiveFlag;
  int m_IsBlackHoleFlag;
  int m_LightScintModelFlag;

  double m_LightBalanceInnerCorr;
  double m_LightBalanceInnerRadius;
  double m_LightBalanceOuterCorr;
  double m_LightBalanceOuterRadius;
};

#endif
