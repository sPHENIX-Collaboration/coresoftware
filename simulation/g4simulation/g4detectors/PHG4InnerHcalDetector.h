#ifndef PHG4InnerHcalDetector_h
#define PHG4InnerHcalDetector_h

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

class PHG4InnerHcalDetector : public PHG4Detector
{
 public:
  typedef CGAL::Exact_circular_kernel_2 Circular_k;
  typedef CGAL::Point_2<Circular_k> Point_2;

  //! constructor
  PHG4InnerHcalDetector(PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4InnerHcalDetector();

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInInnerHcal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  int get_Layer() const { return layer; }
  G4VSolid *ConstructSteelPlate(G4LogicalVolume *hcalenvelope);
  G4VSolid *ConstructScintillatorBox(G4LogicalVolume *hcalenvelope);
  void ShiftSecantToTangent(Point_2 &lowleft, Point_2 &upleft, Point_2 &upright, Point_2 &lowright);

  G4AssemblyVolume *ConstructHcalScintillatorAssembly(G4LogicalVolume *hcalenvelope);
  void ConstructHcalSingleScintillators(G4LogicalVolume *hcalenvelope);
  int CheckTiltAngle() const;
  int ConsistencyCheck() const;
  void SetTiltViaNcross();

 protected:
  int ConstructInnerHcal(G4LogicalVolume *sandwich);
  int DisplayVolume(G4VSolid *volume, G4LogicalVolume *logvol, G4RotationMatrix *rotm = nullptr);
  double x_at_y(Point_2 &p0, Point_2 &p1, double yin);
  PHG4Parameters *params;
  G4AssemblyVolume *scinti_mother_assembly;
  double inner_radius;
  double outer_radius;
  double size_z;
  double scinti_tile_x;
  double scinti_tile_x_lower;
  double scinti_tile_x_upper;
  double scinti_tile_z;
  double scinti_tile_thickness;
  double scinti_inner_gap;
  double scinti_outer_gap;
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
  std::set<G4VPhysicalVolume *> steel_absorber_vec;
  std::vector<G4VSolid *> scinti_tiles_vec;
  std::string scintilogicnameprefix;
};

#endif
