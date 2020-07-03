// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MICROMEGAS_PHG4MICROMEGASSTEPPINGACTION_H
#define G4MICROMEGAS_PHG4MICROMEGASSTEPPINGACTION_H

/*!
 * \file PHG4MicromegasSteppingAction.h
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include <g4main/PHG4SteppingAction.h>

#include <memory>

class PHG4MicromegasDetector;

class G4Step;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHParameters;

class PHG4MicromegasSteppingAction : public PHG4SteppingAction
{

  public:

  //! constructor
  PHG4MicromegasSteppingAction(PHG4MicromegasDetector*, const PHParameters* parameters);

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

  private:

  //! pointer to the detector
  PHG4MicromegasDetector* m_Detector = nullptr;

  const PHParameters* m_Params = nullptr;

  //! pointer to hit container
  PHG4HitContainer* m_hitContainer = nullptr;

  //! running hit
  std::unique_ptr<PHG4Hit> m_hit;

  PHG4HitContainer* m_SaveHitContainer = nullptr;
  G4VPhysicalVolume* m_SaveVolPre = nullptr;
  G4VPhysicalVolume* m_SaveVolPost = nullptr;

  int m_SaveTrackId = -1;
  int m_SavePreStepStatus =  -1;
  int m_SavePostStepStatus =  -1;
  int m_ActiveFlag = 0;
  int m_BlackHoleFlag = 0;
  double m_EdepSum = 0;
  double m_EionSum = 0;
};

#endif // MICROMEGASSTEPPINGACTION_H
