#ifndef PHG4OuterHcalDetector_h
#define PHG4OuterHcalDetector_h

#include "PHG4OuterHcalFieldSetup.h"
#include "PHG4OuterHcalParameters.h"
#include <g4main/PHG4Detector.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4RotationMatrix.hh>

#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>

#include <map>
#include <vector>
#include <set>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;

class PHG4OuterHcalDetector: public PHG4Detector
{
  typedef CGAL::Exact_circular_kernel_2             Circular_k;
  typedef CGAL::Point_2<Circular_k>                 Point_2;

  public:

  //! constructor
  PHG4OuterHcalDetector( PHCompositeNode *Node, PHG4OuterHcalParameters *parameters, const std::string &dnam="BLOCK");

  //! destructor
  virtual ~PHG4OuterHcalDetector();

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  int IsInOuterHcal(G4VPhysicalVolume*) const;
  //@}

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
  int get_Layer() const {return layer;}

  void BlackHole(const int i=1) {blackhole = i;}
  int IsBlackHole() const {return blackhole;}

  protected:
  void AddGeometryNode();
  int ConstructOuterHcal(G4LogicalVolume* sandwich);
  G4VSolid *ConstructHcalSteel(G4LogicalVolume* hcalenvelope);
  G4VSolid *ConstructHcalScintillator(G4LogicalVolume* hcalenvelope);
  int ConstructHcalSingleScintillator(G4LogicalVolume* hcalenvelope);
  G4AssemblyVolume *ConstructHcalScintillatorAssembly(G4LogicalVolume* hcalenvelope);
  int DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix* rotm=NULL);
  PHG4OuterHcalParameters *params;

  // for the initial trapezoid
  G4double steel_rectangle_plate_x; // the rectangle after eta cutout
  G4double steel_plate_x;
  G4double steel_plate_yin;
  G4double steel_plate_yout;
  G4double steel_plate_z;
  G4int n_steel_plates;
  // the scintillator envelop
  G4double scinti_tile_x;
  G4double scinti_tile_y;
  G4double scinti_tile_z;
  G4double scinti_gap;
  G4double scinti_eta_coverage;
  G4double scinti_gap_neighbor;
  G4int n_scinti_tiles;
  G4double gap_between_tiles;
  // the cylinder envelope
  G4double envelope_inner_radius;
  G4double envelope_outer_radius;
  G4double envelope_z;
  //
  G4double tilt_angle;
  // eta at which we cut the steel/scintillator plates to accomodate for magnet related cutout
  G4double etacutline;
  // the box we need to cut out
  G4double cutbox_x;
  G4double cutbox_y;
  G4double cutbox_z;
  G4double cuttrapezoid_x;
  G4double cuttrapezoid_y;
  G4double cuttrapezoid_z_short;
  G4double cuttrapezoid_z_long;
  G4double testbox_x[2];
  G4double testbox_y[2];
  G4double testbox_z[2];
  std::set<G4VPhysicalVolume *>steel_absorber_vec;
  std::set<G4VPhysicalVolume *>scinti_slats_vec;
  /* G4VPhysicalVolume *inner_absorber; */
  /* G4VPhysicalVolume *outer_scinti; */
  /* G4VPhysicalVolume *outer_absorber; */

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
  std::string detector_type;
  std::string superdetector;
  std::vector<G4VSolid *> scinti_tiles_vec; 
  std::string scintilogicnameprefix;
  PHG4OuterHcalFieldSetup * field_setup;
};

#endif
