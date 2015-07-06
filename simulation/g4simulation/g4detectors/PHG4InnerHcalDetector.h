#ifndef PHG4InnerHcalDetector_h
#define PHG4InnerHcalDetector_h

#include "g4main/PHG4Detector.h"

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


class PHG4InnerHcalDetector: public PHG4Detector
{
typedef CGAL::Exact_circular_kernel_2             Circular_k;
typedef CGAL::Point_2<Circular_k>                 Point_2;

  public:

  //! constructor
  PHG4InnerHcalDetector( PHCompositeNode *Node, const std::string &dnam="BLOCK", const int lyr = 0 );

  //! destructor
  virtual ~PHG4InnerHcalDetector(){}

  //! construct
  virtual void Construct( G4LogicalVolume* world );

  //!@name volume accessors
  //@{
  int IsInInnerHcal(G4VPhysicalVolume*) const;
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

  void SetTilt(const double tilt) {tilt_angle = tilt;}

  G4VSolid* ConstructSteelPlate(G4LogicalVolume* hcalenvelope);
  G4VSolid* ConstructScintillatorBox(G4LogicalVolume* hcalenvelope);
  void ShiftSecantToTangent(Point_2 &lowleft, Point_2 &upleft, Point_2 &upright, Point_2 &lowright);

  G4AssemblyVolume *ConstructHcalScintillatorAssembly(G4LogicalVolume* hcalenvelope);
  void ConstructHcalSingleScintillators(G4LogicalVolume* hcalenvelope);
  int CheckTiltAngle() const;

  protected:
  void AddGeometryNode();
  int ConstructInnerHcal(G4LogicalVolume* sandwich);
  int DisplayVolume(G4VSolid *volume,  G4LogicalVolume* logvol, G4RotationMatrix* rotm=NULL);
  G4double x_at_y(Point_2 &p0, Point_2 &p1, G4double yin);
  G4double inner_radius;
  G4double outer_radius;
  G4double size_z;
  G4double scinti_gap;
  G4double tilt_angle;
  int n_scinti_plates;
  G4double scinti_tile_x;
  G4double scinti_tile_y;
  G4double scinti_tile_z;
  int n_scinti_tiles;
  G4double scinti_gap_neighbor;
  G4double scinti_eta_coverage;
  G4double envelope_inner_radius;
  G4double envelope_outer_radius;
  G4double envelope_z;
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
  std::set<G4VPhysicalVolume *>steel_absorber_vec;
  std::vector<G4VSolid *> scinti_tiles_vec; 
  std::string scintilogicnameprefix;
};

#endif
