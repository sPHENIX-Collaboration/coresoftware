// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4BBCSTEPPINGACTION_H
#define G4DETECTORS_PHG4BBCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>  // for string

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4BbcDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class PHG4BbcSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4BbcSteppingAction(PHG4BbcDetector*, const PHParameters*);

  //! destructor
  ~PHG4BbcSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

  void SetHitNodeName(const std::string& type, const std::string& name) override;

 private:
  //! pointer to the detector
  PHG4BbcDetector* m_Detector = nullptr;

  //! pointer to hit container
  PHG4HitContainer* m_HitContainer = nullptr;
  PHG4HitContainer* m_SupportHitContainer = nullptr;
  const PHParameters* m_Params = nullptr;
  PHG4Hit* m_Hit = nullptr;
  PHG4HitContainer* m_SaveHitContainer = nullptr;

  G4VPhysicalVolume* m_SaveVolPre = nullptr;
  G4VPhysicalVolume* m_SaveVolPost = nullptr;

  int m_SaveTrackId = -1;
  int m_SavePreStepStatus = -1;
  int m_SavePostStepStatus = -1;
  int m_ActiveFlag = 0;
  int m_BlackHoleFlag = 0;
  int m_SupportFlag = 0;
  double m_EdepSum = 0.;
  double m_EionSum = 0.;
  double m_PathLen = 0.;

  std::string m_HitNodeName;
  std::string m_SupportNodeName;
};

#endif  // G4DETECTORS_PHG4BBCSTEPPINGACTION_H
