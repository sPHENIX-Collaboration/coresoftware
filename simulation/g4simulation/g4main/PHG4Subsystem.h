// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SUBSYSTEM_H
#define G4MAIN_PHG4SUBSYSTEM_H

#include <fun4all/SubsysReco.h>

#include <iostream>
#include <string>

class PHG4Detector;
class PHG4DisplayAction;
class PHG4EventAction;
class PHG4SteppingAction;
class PHG4TrackingAction;

class PHG4Subsystem: public SubsysReco
{

  public:

  //! constructor
  PHG4Subsystem( const std::string &name = "Generic Subsystem" ): SubsysReco(name),
overlapcheck(false)
  {}

  //! destructor
  virtual ~PHG4Subsystem( void ) {}

  //! event processing
  virtual int process_after_geant(PHCompositeNode *)
  { return 0; }

  //! return pointer to created detector object
  virtual PHG4Detector* GetDetector( void ) const
  { return nullptr; }

  //! return pointer to this subsystem event action
  virtual PHG4EventAction* GetEventAction( void ) const
  { return nullptr; }

  //! return pointer to this subsystem stepping action
  virtual PHG4SteppingAction* GetSteppingAction( void ) const
  { return nullptr; }

  //! return pointer to this subsystem stepping action
  virtual PHG4TrackingAction* GetTrackingAction( void ) const
  { return nullptr; }

  //! return pointer to this subsystem display setting
  virtual PHG4DisplayAction *GetDisplayAction() const
  {return nullptr;}

  void OverlapCheck(const bool chk = true) {overlapcheck = chk;}

  bool CheckOverlap() const {return overlapcheck;}

 private:
  bool overlapcheck;

};

#endif // G4MAIN_PHG4SUBSYSTEM_H
