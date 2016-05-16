// This is the subsystem code (header file) for the 2nd hcal prototype detector
// created on 11/23/2015, HeXC
//
#ifndef PHG4HcalPrototype2Subsystem_h
#define PHG4HcalPrototype2Subsystem_h

#include "g4main/PHG4Subsystem.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4HcalPrototype2Detector;
class PHG4HcalPrototype2SteppingAction;
class PHG4EventAction;

class PHG4HcalPrototype2Subsystem: public PHG4Subsystem
{

  public:

  //! constructor
  PHG4HcalPrototype2Subsystem( const std::string &name = "BLOCK", const int layer = 0 );

  //! destructor
  virtual ~PHG4HcalPrototype2Subsystem( void )
  {}

  //! init
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int Init(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const;

  //void SetXRot(const G4double dbl) {rot_in_x = dbl;}
  void SetYRot(const G4double dbl) {rot_in_y = dbl;}
  void SetZRot(const G4double dbl) {rot_in_z = dbl;}

  void SetMaterial(const std::string &mat) {material = mat;}
  PHG4EventAction* GetEventAction() const {return eventAction_;}
  void SetActive(const int i = 1) {active = i;}
  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() {return superdetector;}

  void BlackHole(const int i=1) {blackhole = i;}

  private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4HcalPrototype2Detector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4HcalPrototype2SteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;
  G4double dimension[3];
  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double rot_in_x;
  G4double rot_in_y;
  G4double rot_in_z;

  G4String material;
  int active;
  int absorberactive;
  int layer;
  int blackhole;
  std::string detector_type;
  std::string superdetector;
};

#endif
