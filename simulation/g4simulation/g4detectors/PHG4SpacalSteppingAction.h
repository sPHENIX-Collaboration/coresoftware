// Tell emacs that this is a C++ source
//  -*- C++ -*-.
/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.1 $$
 * \date $$Date: 2014/03/24 01:36:44 $$
 */

#ifndef G4DETECTORS_PHG4SPACALSTEPPINGACTION_H
#define G4DETECTORS_PHG4SPACALSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4SpacalDetector;
class PHG4Hit;
class PHG4HitContainer;
class PHG4Shower;

class PHG4SpacalSteppingAction : public PHG4SteppingAction
{
 public:
  //! ctor
  explicit PHG4SpacalSteppingAction(PHG4SpacalDetector *);

  //! dtor
  ~PHG4SpacalSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step *, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode *) override;

  double get_zmin() const;

  double get_zmax() const;

  void SetHitNodeName(const std::string &type, const std::string &name) override;

 private:
  //! pointer to the detector
  PHG4SpacalDetector *m_Detector = nullptr;

  //! pointer to hit container
  PHG4HitContainer *m_HitContainer = nullptr;
  PHG4HitContainer *m_AbsorberHitContainer = nullptr;
  PHG4Hit *m_Hit = nullptr;
  PHG4HitContainer *m_CurrentHitContainer = nullptr;
  PHG4Shower *m_CurrentShower = nullptr;
  int m_SaveTrackid = -1;
  int m_SavePostStepStatus = -1;

  std::string m_AbsorberNodeName;
  std::string m_HitNodeName;
};

#endif  // PHG4VHcalSteppingAction_h
