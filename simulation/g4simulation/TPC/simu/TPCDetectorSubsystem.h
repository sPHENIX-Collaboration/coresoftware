#ifndef __TPCDETECTORSUBSYSTEM_H__
#define __TPCDETECTORSUBSYSTEM_H__

#include <string>
#include <g4detectors/PHG4DetectorSubsystem.h>
#include <g4main/PHG4Detector.h>
#include <g4main/PHG4EventAction.h>
#include <g4main/PHG4SteppingAction.h>
#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>
#include "TPCDetector.h"
#include "TPCEventAction.h"
#include "TPCSteppingAction.h"

class TPCDetectorSubsystem: public PHG4DetectorSubsystem {
  public:
  TPCDetectorSubsystem(const std::string &name="", const int lyr=0);
  virtual ~TPCDetectorSubsystem( void ) {}
  void SetDefaultParameters() {}
  int InitRunSubsystem(PHCompositeNode *);
  int process_event(PHCompositeNode *);
  void Print() const;

  PHG4Detector* GetDetector( void ) const {return fDetector;}
  PHG4EventAction* GetEventAction() const {return fEventAction;}
  PHG4SteppingAction* GetSteppingAction( void ) const {return fSteppingAction;}

 private:
  TPCDetector *fDetector;
  TPCEventAction *fEventAction;
  TPCSteppingAction *fSteppingAction;
};

#endif
