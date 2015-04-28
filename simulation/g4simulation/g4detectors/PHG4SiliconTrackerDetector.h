#ifndef PHG4SiliconTrackerDetector_h
#define PHG4SiliconTrackerDetector_h

#include "g4main/PHG4Detector.h"

#include <Geant4/globals.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>

#include <map>
#include <vector>
#include <set>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;


class PHG4SiliconTrackerDetector: public PHG4Detector
{

  public:

  //! constructor
  PHG4SiliconTrackerDetector( PHCompositeNode *Node, const std::string &dnam="BLOCK", const int lyr = 0 );

  //! destructor
  virtual ~PHG4SiliconTrackerDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  int IsInSiliconTracker(G4VPhysicalVolume*) const;
  //@}


  void set_nominal_layer_radius(const G4double radius){layer_nominal_radius = radius * mm;}
  void set_radius_stagger(const G4double stagger){radius_stagger = stagger * mm;}
  void set_N_staggers(const int nstag){N_staggers = nstag;}
  void set_strip_tilt(const G4double tilt){strip_tilt = tilt * rad;}
  void set_option_double_layer(const bool option){option_double_layer = option;}
  void set_N_strips_in_sensor_phi(const int nstrips){N_strips_per_column = nstrips;}
  void set_strip_pitch(const G4double pitch){strip_y = pitch * mm;}
  void set_add_lower_roc(const bool option){add_lower_roc = option;}

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
  int ConstructSiliconTracker(G4LogicalVolume* sandwich);
  int DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix* rotm=NULL);

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
  G4double layer_nominal_radius;
  G4double sensor_z;
  G4double sensor_x;
  G4double strip_x;
  G4double strip_y;
  G4double strip_z;
  G4double roc_y;
  G4double roc_cu_x;
  G4double roc_ins_x;
  G4double tab_z;
  G4double tab_y;
  G4double tab_cu_x;
  G4double tab_ins_x;
  G4double fpga_z;
  G4double fpga_y;
  G4double fpga_x;
  G4double stave_y;
  G4double stave_inner_y;
  G4double stave_x;
  G4double stave_inner_x;
  G4double chip_z;
  G4double chip_y;
  G4double chip_x;
  G4double radius_stagger;
  G4int N_staggers;
  G4double maxrap;
  G4double overlap_fraction;
  G4double strip_tilt;
  G4double sensor_edge_phi;
  G4double sensor_edge_z;
  int N_strips_per_column;
  bool option_double_layer;
  bool add_lower_roc;

  // calculated quantities
  G4double segment_z_spacing;
  G4double segment_phi_spacing;
  G4double strip_z_spacing;
  G4double strip_y_spacing;
  G4double strip_x_offset;
  G4double sensor_y;
  G4double sensor_x_offset;
  G4double sensor_y_offset;
  G4double roc_z;
  G4double stave_z;
  int layer_NPHI;
  int layer_NZ;
  int N_strip_columns;
  int N_sensors_in_layer;

  std::string layer_string;;
  std::string detector_type;
  std::string superdetector;
};

#endif
