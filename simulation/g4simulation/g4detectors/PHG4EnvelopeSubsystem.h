// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4ENVELOPESUBSYSTEM_H
#define G4DETECTORS_PHG4ENVELOPESUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4String.hh>

#include <string>  // for string

class PHG4Detector;
class PHG4EnvelopeDetector;
class PHG4EnvelopeSteppingAction;
class PHG4SteppingAction;
class PHCompositeNode;

class PHG4EnvelopeSubsystem : public PHG4Subsystem
{
 public:
  //Constructor
  PHG4EnvelopeSubsystem(const std::string& name = "ENVELOPE_DEFAULT", const int layer = 0);

  //Destructor
  ~PHG4EnvelopeSubsystem(void) override
  {
  }

  /*
			Creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
			Creates the stepping action and place it on the node tree, under "ACTIONS" node
			Creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
		*/

  int Init(PHCompositeNode*) override;

  //Event Processing
  int process_event(PHCompositeNode*) override;

  //Accessors (reimplemented)
  PHG4Detector* GetDetector(void) const override;
  PHG4SteppingAction* GetSteppingAction(void) const override;

 private:
  //Pointer to Geant4 implementation of detector
  PHG4EnvelopeDetector* detector_;

  //Stepping Action
  PHG4EnvelopeSteppingAction* steppingAction_;

  G4String material;
  int active;

  std::string detector_type;
};

#endif
