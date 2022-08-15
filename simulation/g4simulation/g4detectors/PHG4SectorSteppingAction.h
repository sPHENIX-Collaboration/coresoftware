// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4SECTORSTEPPINGACTION_H
#define G4DETECTORS_PHG4SECTORSTEPPINGACTION_H

#include <g4main/PHG4SteppingAction.h>

class G4Step;
class PHCompositeNode;
class PHG4Hit;
class PHG4HitContainer;
class PHG4SectorDetector;
class PHG4Shower;

class PHG4SectorSteppingAction : public PHG4SteppingAction
{
 public:
  //! constructor
  explicit PHG4SectorSteppingAction(PHG4SectorDetector*);

  //! destructor
  ~PHG4SectorSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  //! pointer to the detector
  PHG4SectorDetector* detector_ = nullptr;

  //! pointer to hit container
  PHG4HitContainer* hits_ = nullptr;
  PHG4Hit* hit = nullptr;
  PHG4Shower* saveshower = nullptr;

  int layer_id = -1;
};

#endif  //__G4PHPHYTHIAREADER_H__
