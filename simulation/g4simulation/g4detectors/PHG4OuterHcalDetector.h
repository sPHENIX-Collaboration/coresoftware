#ifndef PHG4OuterHcalDetector_h
#define PHG4OuterHcalDetector_h

#include "PHG4OuterHcalFieldSetup.h"

#include <g4main/PHG4Detector.h>

// cannot fwd declare G4RotationMatrix, it is a typedef pointing to clhep
#include <Geant4/G4RotationMatrix.hh>

#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>

#include <map>
#include <set>
#include <vector>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHG4Parameters;

class PHG4OuterHcalDetector : public PHG4Detector
{
 public:
  typedef CGAL::Exact_circular_kernel_2 Circular_k;
  typedef CGAL::Point_2<Circular_k> Point_2;
  //! constructor
  PHG4OuterHcalDetector(PHCompositeNode *Node, PHG4Parameters *params, const std::string &dnam = "HCALOUT");

  //! destructor
  virtual ~PHG4OuterHcalDetector();

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInOuterHcal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }
  void ShiftSecantToTangent(Point_2 &lowleft, Point_2 &upleft, Point_2 &upright, Point_2 &lowright);
  int ConsistencyCheck() const;
  void SetTiltViaNcross();
  int CheckTiltAngle() const;
  void ConstructHcalSingleScintillators(G4LogicalVolume *hcalenvelope);
  G4VSolid *ConstructScintillatorBox(G4LogicalVolume *hcalenvelope);

 protected:
  int ConstructOuterHcal(G4LogicalVolume *hcalenvelope);
  G4VSolid *ConstructSteelPlate(G4LogicalVolume *hcalenvelope);
  G4AssemblyVolume *ConstructHcalScintillatorAssembly(G4LogicalVolume *hcalenvelope);
  int DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm = nullptr);
  G4double x_at_y(Point_2 &p0, Point_2 &p1, G4double yin);
  PHG4OuterHcalFieldSetup *field_setup;
  PHG4Parameters *params;
  G4AssemblyVolume *scinti_mother_assembly;
  G4VSolid *steel_cutout_for_magnet;
  double inner_radius;
  double outer_radius;
  double size_z;
  double scinti_tile_x;
  double scinti_tile_x_lower;
  double scinti_tile_x_upper;
  double scinti_tile_z;
  double scinti_tile_thickness;
  double scinti_gap;
  double scinti_inner_radius;
  double scinti_outer_radius;
  double tilt_angle;
  double envelope_inner_radius;
  double envelope_outer_radius;
  double envelope_z;
  double volume_envelope;
  double volume_steel;
  double volume_scintillator;

  int n_scinti_plates;
  int n_scinti_tiles;

  int active;
  int absorberactive;

  int layer;
  std::string detector_type;
  std::string superdetector;
  std::string scintilogicnameprefix;
  std::vector<G4VSolid *> scinti_tiles_vec;
  std::set<G4VPhysicalVolume *> steel_absorber_vec;
};

#endif
