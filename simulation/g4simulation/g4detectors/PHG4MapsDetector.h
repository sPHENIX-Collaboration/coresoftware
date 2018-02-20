#ifndef PHG4MapsDetector_h
#define PHG4MapsDetector_h

#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>

#include <map>
#include <vector>
#include <set>
#include <string>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHParameters;

class PHG4MapsDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4MapsDetector( PHCompositeNode *Node,  PHParameters *parameters, const std::string &dnam="BLOCK" );

  //! destructor
  virtual ~PHG4MapsDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  int IsInMaps(G4VPhysicalVolume*) const;
  //@}

  void set_stave_type(const int st){stave_type = st;}
  void set_nominal_layer_radius(const G4double radius){layer_nominal_radius = radius * mm;}
  void set_pixel_x(const double pixel_x_in) {pixel_x = pixel_x_in;}
  void set_pixel_z(const double pixel_z_in) {pixel_z = pixel_z_in;}
  void set_pixel_thickness(const double pixel_thickness_in) {pixel_thickness = pixel_thickness_in;}
  void SetPlaceZ(const G4double place_z) {place_in_z = place_z*cm;}
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z)
  {
    place_in_x = place_x*cm;
    place_in_y = place_y*cm;
    place_in_z = place_z*cm;
  }
  void SetXRot(const G4double angle) {x_rot = angle*rad;}
  void SetYRot(const G4double angle) {y_rot = angle*rad;}
  void SetZRot(const G4double angle) {z_rot = angle*rad;}
  void SetActive(const int i = 1) {active = i;}
  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  int IsActive() const {return active;}
 void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() const {return superdetector;}
  void Detector(const std::string &name) {detector_type = name;}
  const std::string Detector() const {return detector_type;}
  int get_Layer() const {return layer;}

  void BlackHole(const int i=1) {blackhole = i;}
  int IsBlackHole() const {return blackhole;}

  private:
  void AddGeometryNode();
  int ConstructMaps(G4LogicalVolume* sandwich);
  void SetDisplayProperty( G4AssemblyVolume* av);
  void SetDisplayProperty( G4LogicalVolume* lv);

  // the cylinder envelope
  G4double envelope_inner_radius;
  G4double envelope_outer_radius;
  G4double envelope_z;
  //
   G4double place_in_x;
   G4double place_in_y;
   G4double place_in_z;
   G4double x_rot;
   G4double y_rot;
   G4double z_rot;
  int active;
  int absorberactive;
  int layer;
  int blackhole;

  // setup parameters
  int stave_type;
  G4double layer_nominal_radius;
  int N_staves;
  G4double phistep;
  G4double phitilt;
  double pixel_x;
  double pixel_z;
  double pixel_thickness;

  // calculated quantities

  std::string layer_string;;
  std::string detector_type;
  std::string superdetector;
  std::string stave_geometry_file;

//  int verbosity;
};

#endif
