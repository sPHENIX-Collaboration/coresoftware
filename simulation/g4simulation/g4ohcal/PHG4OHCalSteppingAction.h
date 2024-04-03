// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4OHCAL_PHG4OHCALSTEPPINGACTION_H
#define G4OHCAL_PHG4OHCALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <string>  // for string

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class TowerInfoContainer;
class PHG4OHCalDetector;
class PHParameters;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class TH2;

class PHG4OHCalSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  PHG4OHCalSteppingAction(PHG4OHCalDetector *, PHParameters *parameters);

  //! destructor
  ~PHG4OHCalSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  int InitWithNode(PHCompositeNode *topNode) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

  void FieldChecker(const G4Step *);
  void EnableFieldChecker(const int i = 1) { m_EnableFieldCheckerFlag = i; }
  void CreateNodeTree(PHCompositeNode *topNode);

 private:
  bool NoHitSteppingAction(const G4Step *aStep);
  //! pointer to the detector
  PHG4OHCalDetector *m_Detector{nullptr};

  //! efficiency maps from Mephi
  TH2 *m_MapCorrHist[24]{0};
  TH2 *m_MapCorrHistChim[24]{0};

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
  int m_EnableFieldCheckerFlag{-1};

  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActiveFlag{0};
  int m_IsBlackHoleFlag{0};
  int m_NScintiPlates{-1};
  int m_LightScintModelFlag{0};
  bool m_doG4Hit{true};
  double m_tmin{-20.};
  double m_tmax{60.};
  double m_dt{100.};
  std::string m_AbsorberNodeName;
  std::string m_HitNodeName;

  TowerInfoContainer *m_CaloInfoContainer{nullptr};
};

#endif  // G4OHCAL_PHG4OHCALSTEPPINGACTION_H
