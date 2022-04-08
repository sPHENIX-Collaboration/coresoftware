// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTSTEPPINGACTION_H
#define G4INTT_PHG4INTTSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <map>
#include <string>
#include <utility>  // for pair
#include <vector>

class G4Step;
class PHCompositeNode;
class PHG4InttDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;
class PHParametersContainer;

class PHG4InttSteppingAction : public PHG4SteppingAction
{
 public:
  PHG4InttSteppingAction(PHG4InttDetector *, const PHParametersContainer *parameters, const std::pair<std::vector<std::pair<int, int>>::const_iterator, std::vector<std::pair<int, int>>::const_iterator> &layer_begin_end);

  ~PHG4InttSteppingAction() override;

  bool UserSteppingAction(const G4Step *, bool) override;

  void SetInterfacePointers(PHCompositeNode *) override;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

 private:
  //! pointer to the detector
  PHG4InttDetector *m_Detector = nullptr;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer = nullptr;
  PHG4HitContainer *m_AbsorberHitContainer = nullptr;
  PHG4Hit *m_Hit = nullptr;
  PHG4HitContainer *m_SaveHitContainer = nullptr;
  PHG4Shower *m_SaveShower = nullptr;
  const PHParametersContainer *m_ParamsContainer = nullptr;

  std::map<int, int> m_InttToTrackerLayerMap;
  std::map<int, int> m_LadderTypeMap;
  std::map<int, double> m_StripYMap;
  std::map<int, std::pair<double, double>> m_StripZMap;
  std::map<int, int> m_nStripsPhiCell;
  std::map<int, std::pair<int, int>> m_nStripsZSensor;

  std::map<int, int> m_IsActiveMap;
  std::map<int, int> m_IsBlackHoleMap;
  std::string m_AbsorberNodeName;
  std::string m_HitNodeName;
};

#endif  // G4INTT_PHG4INTTSTEPPINGACTION_H
