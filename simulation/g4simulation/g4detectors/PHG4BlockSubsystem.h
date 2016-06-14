#ifndef PHG4BlockSubsystem_h
#define PHG4BlockSubsystem_h

#include "PHG4DetectorSubsystem.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

class PHG4BlockDetector;
class PHG4BlockSteppingAction;
class PHG4EventAction;

class PHG4BlockSubsystem: public PHG4DetectorSubsystem
{
 public:

  //! constructor
  PHG4BlockSubsystem( const std::string &name = "BLOCK", const int layer = 0 );

  //! destructor
  virtual ~PHG4BlockSubsystem( void )
  {}

  //! InitRunSubsystem
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode *);

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
  void SetZRot(const G4double d) {_rot_in_z = d;}
  PHG4EventAction* GetEventAction() const {return _eventAction;}
  void SuperDetector(const std::string &name) {_superdetector = name;}
  const std::string SuperDetector() {return _superdetector;}

  void BlackHole(const int i=1) {_blackhole = i;}
  void UseG4Steps(const int i = 1);
  void UseIonizationEnergy(const int i = 1);

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4BlockDetector* _detector;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SteppingAction* _steppingAction;
  PHG4EventAction *_eventAction;
  G4double _dimension[3];
  G4double _rot_in_z;

  int _layer;
  int _blackhole;
  int _use_g4_steps;
  int _use_ionisation_energy;
  std::string _detector_type;
  std::string _superdetector;
};

#endif
