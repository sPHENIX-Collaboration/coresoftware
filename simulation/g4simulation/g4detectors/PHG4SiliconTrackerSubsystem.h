#ifndef PHG4SiliconTrackerSubsystem_h
#define PHG4SiliconTrackerSubsystem_h

#include "g4main/PHG4Subsystem.h"

#include <Geant4/G4Types.hh>
#include <Geant4/G4String.hh>

#include <utility>
#include <vector>

class PHG4SiliconTrackerDetector;
class PHG4SiliconTrackerSteppingAction;
class PHG4EventAction;

typedef std::vector<std::pair<int, int>> vpair;

class PHG4SiliconTrackerSubsystem: public PHG4Subsystem
{
  public:

  //! constructor
  PHG4SiliconTrackerSubsystem(const std::string &name = "BLOCK", const vpair &layerconfig=vpair(0));

  //! destructor
  virtual ~PHG4SiliconTrackerSubsystem( void )
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

  // These are passed on to the detector class, which expects mm and radians
  void set_nominal_layer_radius(const G4double radius){layer_nominal_radius = radius;}
  void set_radius_stagger(const G4double stagger){radius_stagger = stagger;}
  void set_N_staggers(const int nstag){N_staggers = nstag;}
  void set_N_strips_in_sensor_phi(const int nstrips){N_strips_per_column = nstrips;}
  void set_strip_tilt(const G4double tilt){strip_tilt = tilt;}
  void set_option_double_layer(const bool option){option_double_layer = option;}
  void set_add_lower_roc(const bool option){add_lower_roc = option;}

  void SetSize(const G4double sizex, const G4double sizey, const G4double sizez)
     {dimension[0] = sizex; dimension[1] = sizey; dimension[2] = sizez;}
  void SetPlaceZ(const G4double dbl) {place_in_z = dbl;}
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
    place_in_x = place_x;
    place_in_y = place_y;
    place_in_z = place_z;
  }
  void SetXRot(const G4double dbl) {rot_in_x = dbl;}
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
  PHG4SiliconTrackerDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4SiliconTrackerSteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;
  G4double dimension[3];
  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double rot_in_x;
  G4double rot_in_y;
  G4double rot_in_z;

  // These are passed on to the detector class
  G4double layer_nominal_radius;
  G4double radius_stagger;
  G4int N_staggers;
  G4double sensor_y;
  G4double strip_tilt;
  bool option_double_layer;
  bool add_lower_roc;
  int N_strips_per_column;

  G4String material;
  int active;
  int absorberactive;
  std::vector< std::pair<int,int> > layerconfig_;
  int blackhole;
  std::string detector_type;
  std::string superdetector;
};

#endif
