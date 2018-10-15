// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROTOTYPE2OUTERHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4PROTOTYPE2OUTERHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class PHG4Prototype2OuterHcalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4Prototype2OuterHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4Prototype2OuterHcalSteppingAction(PHG4Prototype2OuterHcalDetector *, const PHParameters *parameters);

  //! dtor
  virtual ~PHG4Prototype2OuterHcalSteppingAction();

  //! stepping action
  virtual bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  virtual void SetInterfacePointers(PHCompositeNode *);

  double GetLightCorrection(const double r) const;

 private:
  //! pointer to the detector
  PHG4Prototype2OuterHcalDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *m_Hits;
  PHG4HitContainer *m_AbsorberHits;
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
