// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4CONESUBSYSTEM_H
#define G4DETECTORS_PHG4CONESUBSYSTEM_H

#include <g4main/PHG4Subsystem.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <string>                  // for string

class PHCompositeNode;
class PHG4ConeDetector;
class PHG4ConeSteppingAction;
class PHG4Detector;
class PHG4EventAction;
class PHG4SteppingAction;

class PHG4ConeSubsystem: public PHG4Subsystem
{

  public:

  //! constructor
  PHG4ConeSubsystem( const std::string &name = "CONE", const int layer = 0 );

  //! destructor
  virtual ~PHG4ConeSubsystem( void )
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

  //!set inner and outter radius1
  void SetR1(const G4double min, const G4double max)
  {rMin1 = min; rMax1=max;}

  //!set inner and outter radius2
  void SetR2(const G4double min, const G4double max)
  {rMin2 = min; rMax2=max;}

  //! set length in Z
  void SetZlength(const G4double a)
  {dZ = a;}

  //! set phi offset and extention
  void SetPhi(const G4double a, const G4double b)
  {sPhi = a; dPhi = b;}

  //! set rmaximum and minimums according to the eta range 
  void Set_eta_range(G4double etaMin, G4double etaMax);

  void SetPlaceZ(const G4double dbl) {place_in_z = dbl;}
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
    place_in_x = place_x;
    place_in_y = place_y;
    place_in_z = place_z;
  }
  void SetZRot(const G4double dbl) {rot_in_z = dbl;}
  void SetMaterial(const std::string &mat) {material = mat;}
  PHG4EventAction* GetEventAction() const {return eventAction_;}
  void SetActive(const int i = 1) {active = i;}
  void SuperDetector(const std::string &name) {superdetector = name;}

  private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4ConeDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4ConeSteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;

  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double rot_in_z;
  G4double rMin1;
  G4double rMax1;
  G4double rMin2;
  G4double rMax2;
  G4double dZ;
  G4double sPhi;
  G4double dPhi;
  G4String material;
  int active;
  int layer;
  std::string detector_type;
  std::string superdetector;
};

#endif
