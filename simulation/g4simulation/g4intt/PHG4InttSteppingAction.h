// Tell emacs that this is a C++ source
// This file is really -*- C++ -*-.
#ifndef G4INTT_PHG4INTTSTEPPINGACTION_H
#define G4INTT_PHG4INTTSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

#include <map>
#include <vector>
#include <utility>                      // for pair

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

  virtual ~PHG4InttSteppingAction();

  virtual bool UserSteppingAction(const G4Step *, bool);

  virtual void SetInterfacePointers(PHCompositeNode *);

 private:
  //! pointer to the detector
  PHG4InttDetector *m_Detector;

  //! pointer to hit container
  PHG4HitContainer *m_Hits;
  PHG4HitContainer *m_AbsorberHits;
  PHG4Hit *m_Hit;
  PHG4HitContainer *m_SaveHitContainer;
  PHG4Shower *m_SaveShower;
  const PHParametersContainer *m_ParamsContainer;

  std::map<int, int> m_InttToTrackerLayerMap;
  std::map<int, int> m_LadderTypeMap;
  std::map<int, double> m_StripYMap;
  std::map<int, std::pair<double, double>> m_StripZMap;
  std::map<int, int> m_nStripsPhiCell;
  std::map<int, std::pair<int, int>> m_nStripsZSensor;

  std::map<int, int> m_IsActiveMap;
  std::map<int, int> m_IsBlackHoleMap;
};

#endif  // G4INTT_PHG4INTTSTEPPINGACTION_H
