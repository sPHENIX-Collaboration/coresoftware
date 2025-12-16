// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4OUTERHCALDETECTOR_H
#define G4DETECTORS_PHG4OUTERHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>  // for G4double

#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>

#include <limits>
#include <map>
#include <set>
#include <string>   // for string
#include <utility>  // for pair
#include <vector>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHCompositeNode;
class PHG4OuterHcalDisplayAction;
class PHG4OuterHcalFieldSetup;
class PHParameters;
class PHG4Subsystem;
class RawTowerGeomContainer;

class PHG4OuterHcalDetector : public PHG4Detector
{
 public:
  typedef CGAL::Exact_circular_kernel_2 Circular_k;
  typedef CGAL::Point_2<Circular_k> Point_2;
  //! constructor
  PHG4OuterHcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam);

  //! destructor
  ~PHG4OuterHcalDetector() override;

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInOuterHcal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string &SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  void ShiftSecantToTangent(Point_2 &lowleft, Point_2 &upleft, Point_2 &upright, Point_2 &lowright) const;
  int ConsistencyCheck() const;
  void SetTiltViaNcross();
  int CheckTiltAngle() const;
  void ConstructHcalSingleScintillators(G4LogicalVolume *hcalenvelope);
  G4VSolid *ConstructScintillatorBox(G4LogicalVolume *hcalenvelope);
  std::pair<int, int> GetLayerTowerId(G4VPhysicalVolume *volume) const;

 protected:
  void AddGeometryNode();
  int ConstructOuterHcal(G4LogicalVolume *hcalenvelope);
  G4VSolid *ConstructSteelPlate(G4LogicalVolume *hcalenvelope);
  G4AssemblyVolume *ConstructHcalScintillatorAssembly(G4LogicalVolume *hcalenvelope);
  G4double x_at_y(Point_2 &p0, Point_2 &p1, G4double yin);
  PHG4OuterHcalDisplayAction *m_DisplayAction {nullptr};
  PHG4OuterHcalFieldSetup *m_FieldSetup {nullptr};
  PHParameters *m_Params {nullptr};
  G4AssemblyVolume *m_ScintiMotherAssembly {nullptr};
  G4VSolid *m_SteelCutoutForMagnetG4Solid {nullptr};
  double m_InnerRadius {std::numeric_limits<double>::quiet_NaN()};
  double m_OuterRadius {std::numeric_limits<double>::quiet_NaN()};
  double m_SizeZ {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiTileX {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiTileXLower {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiTileXUpper {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiTileZ {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiTileThickness {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiGap {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiInnerRadius {std::numeric_limits<double>::quiet_NaN()};
  double m_ScintiOuterRadius {std::numeric_limits<double>::quiet_NaN()};
  double m_TiltAngle {std::numeric_limits<double>::quiet_NaN()};
  double m_EnvelopeInnerRadius {std::numeric_limits<double>::quiet_NaN()};
  double m_EnvelopeOuterRadius {std::numeric_limits<double>::quiet_NaN()};
  double m_EnvelopeZ {std::numeric_limits<double>::quiet_NaN()};
  double m_VolumeEnvelope {std::numeric_limits<double>::quiet_NaN()};
  double m_VolumeSteel {std::numeric_limits<double>::quiet_NaN()};
  double m_VolumeScintillator {std::numeric_limits<double>::quiet_NaN()};

  int m_NumScintiPlates {0};
  int m_NumScintiTiles {0};

  int m_ActiveFlag {0};
  int m_AbsorberActiveFlag {0};

  int m_Layer {0};
  std::string m_SuperDetector;
  std::string m_ScintiLogicNamePrefix;
  std::vector<G4VSolid *> m_ScintiTilesVec;
  std::set<G4VPhysicalVolume *> m_SteelAbsorberVec;
  std::map<G4VPhysicalVolume *, std::pair<int, int>> m_ScintiTilePhysVolMap;
  RawTowerGeomContainer *m_RawTowerGeom {nullptr};
  std::string m_TowerGeomNodeName;
};

#endif
