#ifndef G4TPC_PHG4VTPCSTEPPINGACTION_H
#define G4TPC_PHG4VTPCSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4TpcDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParameters;

class PHG4TpcSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4TpcSteppingAction(PHG4TpcDetector *, const PHParameters *parameters);

  //! destructor
  ~PHG4TpcSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

 private:
  //! pointer to the detector
  PHG4TpcDetector *m_Detector = nullptr;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer = nullptr;
  PHG4HitContainer *m_AbsorberHitContainer = nullptr;
  PHG4Hit *m_Hit = nullptr;
  const PHParameters *m_Params = nullptr;
  PHG4HitContainer *m_CurrentHitContainer = nullptr;
  PHG4Shower *m_Shower = nullptr;
  G4VPhysicalVolume *m_SaveVolPre = nullptr;
  G4VPhysicalVolume *m_SaveVolPost = nullptr;
  int m_SaveTrackId = -1;
  int m_SavePreStepStatus = -1;
  int m_SavePostStepStatus = -1;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsBlackHoleFlag = 0;
  int m_UseG4StepsFlag = 0;

  std::string m_HitNodeName;
  std::string m_AbsorberNodeName;
};

#endif  // G4TPC_PHG4TPCSTEPPINGACTION_H
