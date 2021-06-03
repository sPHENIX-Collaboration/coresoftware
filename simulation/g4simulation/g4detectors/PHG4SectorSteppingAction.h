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
  PHG4SectorSteppingAction(PHG4SectorDetector*);

  //! destructor
  ~PHG4SectorSteppingAction() override;

  //! stepping action
  bool UserSteppingAction(const G4Step*, bool) override;

  //! reimplemented from base class
  void SetInterfacePointers(PHCompositeNode*) override;

 private:
  //! pointer to the detector
  PHG4SectorDetector* detector_;

  //! pointer to hit container
  PHG4HitContainer* hits_;
  PHG4Hit* hit;
  PHG4Shower* saveshower;

  int layer_id;
};

#endif  //__G4PHPHYTHIAREADER_H__
