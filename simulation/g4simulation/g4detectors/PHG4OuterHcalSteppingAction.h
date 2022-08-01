// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4OUTERHCALSTEPPINGACTION_H
#define G4DETECTORS_PHG4OUTERHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4OuterHcalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class TH2;

class PHG4OuterHcalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4OuterHcalSteppingAction(PHG4OuterHcalDetector *, const PHParameters *parameters);

  //! destructor
  ~PHG4OuterHcalSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  int Init() override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  void FieldChecker(const G4Step *);
  void EnableFieldChecker(const int i = 1) { m_EnableFieldCheckerFlag = i; }

 private:
  //! pointer to the detector
  PHG4OuterHcalDetector *m_Detector = nullptr;

  //! efficiency maps from Mephi
  TH2 *m_MapCorrHist = nullptr;

  //! pointer to hit container
  PHG4HitContainer *m_Hits = nullptr;
  PHG4HitContainer *m_AbsorberHits = nullptr;
  PHG4Hit *m_Hit = nullptr;
  const PHParameters *m_Params = nullptr;
  PHG4HitContainer *m_SaveHitContainer = nullptr;
  PHG4Shower *m_SaveShower = nullptr;
  G4VPhysicalVolume *m_SaveVolPre = nullptr;
  G4VPhysicalVolume *m_SaveVolPost = nullptr;
  int m_SaveTrackId = -1;
  int m_SavePreStepStatus = -1;
  int m_SavePostStepStatus = -1;
  int m_EnableFieldCheckerFlag = -1;

  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActiveFlag = -1;
  int m_IsBlackHoleFlag = -1;
  int m_NScintiPlates = -1;
  int m_LightScintModelFlag = 0;
};

#endif  // G4DETECTORS_PHG4OUTERHCALSTEPPINGACTION_H
