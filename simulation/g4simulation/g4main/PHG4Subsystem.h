// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4SUBSYSTEM_H
#define G4MAIN_PHG4SUBSYSTEM_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <iostream>
#include <string>

class G4LogicalVolume;
class PHG4Detector;
class PHG4DisplayAction;
class PHG4EventAction;
class PHG4SteppingAction;
class PHG4TrackingAction;

class PHG4Subsystem: public SubsysReco
{

  public:

  //! constructor
  PHG4Subsystem( const std::string &name = "Generic Subsystem" ): 
SubsysReco(name),
  m_MyMotherSubsystem(nullptr),
m_MyLogicalVolume(nullptr),
overlapcheck(false)
  {}

  //! destructor
  virtual ~PHG4Subsystem( void ) {}

  //! event processing
  virtual int process_after_geant(PHCompositeNode *)
  { return Fun4AllReturnCodes::EVENT_OK; }

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

  void SetMotherSubsystem(PHG4Subsystem *subsys) {m_MyMotherSubsystem = subsys;}
  PHG4Subsystem *GetMotherSubsystem() const {return m_MyMotherSubsystem;}

  void SetLogicalVolume(G4LogicalVolume *vol) {m_MyLogicalVolume = vol;}
  G4LogicalVolume *GetLogicalVolume() const {return m_MyLogicalVolume;}

 private:
  PHG4Subsystem *m_MyMotherSubsystem;
  G4LogicalVolume *m_MyLogicalVolume;
  bool overlapcheck;

};

#endif // G4MAIN_PHG4SUBSYSTEM_H
