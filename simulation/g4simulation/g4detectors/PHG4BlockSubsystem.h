#ifndef PHG4BlockSubsystem_h
#define PHG4BlockSubsystem_h

#include "g4main/PHG4Subsystem.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4BlockDetector;
class PHG4BlockSteppingAction;
class PHG4EventAction;

class PHG4BlockSubsystem: public PHG4Subsystem
{
 public:

  //! constructor
  PHG4BlockSubsystem( const std::string &name = "BLOCK", const int layer = 0 );

  //! destructor
  virtual ~PHG4BlockSubsystem( void )
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

  void SetSize(const G4double sizex, const G4double sizey, const G4double sizez)
    {_dimension[0] = sizex; _dimension[1] = sizey; _dimension[2] = sizez;}
  void SetCenterZ(const G4double dbl) {_center_in_z = dbl;}
  void SetCenter(const G4double center_x, const G4double center_y, const G4double center_z)
  {
    _center_in_x = center_x;
    _center_in_y = center_y;
    _center_in_z = center_z;
  }

  void SetZRot(const G4double d) {_rot_in_z = d;}
  void SetMaterial(const std::string &mat) {_material = mat;}
  PHG4EventAction* GetEventAction() const {return _eventAction;}
  void SetActive(const int i = 1) {_active = i;}
  void SuperDetector(const std::string &name) {_superdetector = name;}
  const std::string SuperDetector() {return _superdetector;}

  void BlackHole(const int i=1) {_blackhole = i;}
  void UseG4Steps(const int i = 1);
  void UseIonizationEnergy(const int i = 1);

 private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4BlockDetector* _detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4BlockSteppingAction* _steppingAction;
  PHG4EventAction *_eventAction;
  G4double _dimension[3];
  G4double _center_in_x;
  G4double _center_in_y;
  G4double _center_in_z;
  G4double _rot_in_z;

  G4String _material;
  int _active;
  int _layer;
  int _blackhole;
  int _use_g4_steps;
  int _use_ionisation_energy;
  std::string _detector_type;
  std::string _superdetector;
};

#endif
