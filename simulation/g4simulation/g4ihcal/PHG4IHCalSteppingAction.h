// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4IHCAL_PHG4IHCALSTEPPINGACTION_H
#define G4IHCAL_PHG4IHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>  // for string

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4IHCalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class TH2;

class PHG4IHCalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4IHCalSteppingAction(PHG4IHCalDetector *, const PHParameters *parameters);

  //! destructor
  ~PHG4IHCalSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  int Init() override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

 private:
  //! pointer to the detector
  PHG4IHCalDetector *m_Detector = nullptr;

  //! efficiency maps from Mephi
  TH2 *m_MapCorrHist = nullptr;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer = nullptr;
  PHG4HitContainer *m_AbsorberHitContainer = nullptr;
  PHG4Hit *m_Hit = nullptr;
  const PHParameters *m_Params = nullptr;
  PHG4HitContainer *m_SaveHitContainer = nullptr;
  PHG4Shower *m_SaveShower = nullptr;
  G4VPhysicalVolume *m_SaveVolPre = nullptr;
  G4VPhysicalVolume *m_SaveVolPost = nullptr;
  int m_SaveTrackId = -1;
  int m_SavePreStepStatus = -1;
  int m_SavePostStepStatus = -1;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActive = 0;
  int m_IsBlackHole = 0;
  int m_LightScintModelFlag = 0;
  std::string m_AbsorberNodeName;
  std::string m_HitNodeName;
};

#endif  // G4IHCAL_PHG4IHCALSTEPPINGACTION_H
