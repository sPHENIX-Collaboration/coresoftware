// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROTOTYPE3INNERHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4PROTOTYPE3INNERHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Prototype3InnerHcalDetector;
class PHG4Shower;
class PHParameters;

class PHG4Prototype3InnerHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4Prototype3InnerHcalSteppingAction(PHG4Prototype3InnerHcalDetector *, const PHParameters *parameters);

  //! dtor
  virtual ~PHG4Prototype3InnerHcalSteppingAction();

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool);

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4Prototype3InnerHcalDetector *m_Detector;

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
  int m_AbsorberTruth;
  int m_IsActive;
  int m_IsBlackHole;
  int m_LightScintModel;

};

#endif  // G4DETECTORS_PHG4PROTOTYPE2INNERHCALSTEPPINGACTION_H
