// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHG4TPCENDCAPSTEPPINGACTION_H
#define PHG4TPCENDCAPSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>  // for string

class PHG4TpcEndCapDetector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class PHG4TpcEndCapSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4TpcEndCapSteppingAction(PHG4TpcEndCapDetector*, const PHParameters* parameters);

  //! destructor
  ~PHG4TpcEndCapSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

  void SetHitNodeName(const std::string& type, const std::string& name) override;

 private:
  //! pointer to the detector
  PHG4TpcEndCapDetector* m_Detector = nullptr;
  const PHParameters* m_Params = nullptr;
  //! pointer to hit container
  PHG4HitContainer* m_HitContainer = nullptr;
  PHG4Hit* m_Hit = nullptr;
  PHG4HitContainer* m_SaveHitContainer = nullptr;
  G4VPhysicalVolume* m_SaveVolPre = nullptr;
  G4VPhysicalVolume* m_SaveVolPost = nullptr;

  int m_SaveTrackId = -1;
  int m_SavePreStepStatus = -1;
  int m_SavePostStepStatus = -1;
  int m_ActiveFlag = 0;
  int m_BlackHoleFlag = 0;
  double m_EdepSum = 0.;
  double m_EionSum = 0.;

  std::string m_HitNodeName;
};

#endif  // PHG4TPCENDCAPSTEPPINGACTION_H
