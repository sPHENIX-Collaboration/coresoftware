// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4IHCAL_PHG4IHCALSTEPPINGACTION_H
#define G4IHCAL_PHG4IHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>  // for string

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class TowerInfoContainer;
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
  PHG4IHCalSteppingAction(PHG4IHCalDetector *, PHParameters *parameters);

  //! destructor
  ~PHG4IHCalSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  int InitWithNode(PHCompositeNode *topNode) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

  void CreateNodeTree(PHCompositeNode *topNode);

 private:
  bool NoHitSteppingAction(const G4Step *aStep);
  //! pointer to the detector
  PHG4IHCalDetector *m_Detector{nullptr};

  //! efficiency maps from Mephi
  TH2 *m_MapCorrHist{nullptr};

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer{nullptr};
  PHG4HitContainer *m_AbsorberHitContainer{nullptr};
  PHG4Hit *m_Hit{nullptr};
  PHParameters *m_Params{nullptr};
  PHG4HitContainer *m_SaveHitContainer{nullptr};
  PHG4Shower *m_SaveShower{nullptr};
  G4VPhysicalVolume *m_SaveVolPre{nullptr};
  G4VPhysicalVolume *m_SaveVolPost{nullptr};
  int m_SaveTrackId{-1};
  int m_SavePreStepStatus{-1};
  int m_SavePostStepStatus{-1};
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActive{0};
  int m_IsBlackHole{0};
  int m_LightScintModelFlag{0};
  bool m_doG4Hit{true};
  double m_tmin{-20.};
  double m_tmax{60.};
  double m_dt{100.};

  std::string m_AbsorberNodeName;
  std::string m_HitNodeName;

  TowerInfoContainer *m_CaloInfoContainer{nullptr};
};

#endif  // G4IHCAL_PHG4IHCALSTEPPINGACTION_H
