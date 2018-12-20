#ifndef PHG4MVTXSubsystem_h
#define PHG4MVTXSubsystem_h

#include <g4detectors/PHG4DetectorSubsystem.h>

#include <Geant4/G4String.hh>
#include <Geant4/G4Types.hh>

class PHG4MVTXDetector;
class PHG4MVTXSteppingAction;
class PHG4EventAction;

class PHG4MVTXSubsystem : public PHG4DetectorSubsystem
{
 public:
  //! constructor
  PHG4MVTXSubsystem(const std::string& name = "BLOCK", const int layer = 0, int stave_type = 0);

  //! destructor
  virtual ~PHG4MVTXSubsystem(void)
  {
  }

  //! init
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int InitRunSubsystem(PHCompositeNode*);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode*);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector(void) const;
  virtual PHG4SteppingAction* GetSteppingAction(void) const;

  // These are passed on to the detector class, which expects mm and radians
  //  void set_nominal_layer_radius(const G4double radius){layer_nominal_radius = radius;}
  //  void set_stave_type(const int stype){stave_type = stype;}
  //  void set_pixel_x(const double pixel_x_in) {pixel_x = pixel_x_in;}
  //  void set_pixel_z(const double pixel_z_in) {pixel_z = pixel_z_in;}
  //  void set_pixel_thickness(const double pixel_thickness_in) {pixel_thickness = pixel_thickness_in;}
  //  void SetSize(const G4double sizex, const G4double sizey, const G4double sizez)
  //     {dimension[0] = sizex; dimension[1] = sizey; dimension[2] = sizez;}
  //  void SetPlaceZ(const G4double dbl) {place_in_z = dbl;}
  //  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  //  {
  //    place_in_x = place_x;
  //    place_in_y = place_y;
  //    place_in_z = place_z;
  //  }
  //  void SetXRot(const G4double dbl) {rot_in_x = dbl;}
  //  void SetYRot(const G4double dbl) {rot_in_y = dbl;}
  //  void SetZRot(const G4double dbl) {rot_in_z = dbl;}
  //  void SetMaterial(const std::string &mat) {material = mat;}
  PHG4EventAction* GetEventAction() const { return eventAction_; }
  //  void SetActive(const int i = 1) {active = i;}
  //  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  void SuperDetector(const std::string& name) { superdetector = name; }
  const std::string SuperDetector() { return superdetector; }

  //  void BlackHole(const int i=1) {blackhole = i;}

 private:
  void SetDefaultParameters();

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4MVTXDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4MVTXSteppingAction* steppingAction_;
  PHG4EventAction* eventAction_;
  //  G4double dimension[3];
  //  G4double place_in_x;
  //  G4double place_in_y;
  //  G4double place_in_z;
  //  G4double rot_in_x;
  //  G4double rot_in_y;
  //  G4double rot_in_z;

  // These are passed on to the detector class
  G4double layer_nominal_radius;
  int layer;
  int stave_type;

  //  double pixel_x;
  //  double pixel_z;
  //  double pixel_thickness;
  //
  //  G4String material;
  //  int active;
  //  int absorberactive;
  //  int blackhole;
  std::string detector_type;
  std::string superdetector;
};

#endif
