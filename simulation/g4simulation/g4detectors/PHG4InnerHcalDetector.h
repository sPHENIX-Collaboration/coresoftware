// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4INNERHCALDETECTOR_H
#define G4DETECTORS_PHG4INNERHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <CGAL/Cartesian.h>  // for Cartesian_base_ref_count::...
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Point_2.h>  // for Point_2
#pragma GCC diagnostic pop

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
class PHG4InnerHcalDisplayAction;
class PHParameters;
class PHG4Subsystem;

class PHG4InnerHcalDetector : public PHG4Detector
{
 public:
  typedef CGAL::Exact_circular_kernel_2 Circular_k;
  typedef CGAL::Point_2<Circular_k> Point_2;

  //! constructor
  PHG4InnerHcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  ~PHG4InnerHcalDetector() override;

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInInnerHcal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  G4VSolid *ConstructSteelPlate(G4LogicalVolume *hcalenvelope);
  G4VSolid *ConstructScintillatorBox(G4LogicalVolume *hcalenvelope);
  void ShiftSecantToTangent(Point_2 &lowleft, Point_2 &upleft, Point_2 &upright, Point_2 &lowright);

  G4AssemblyVolume *ConstructHcalScintillatorAssembly(G4LogicalVolume *hcalenvelope);
  void ConstructHcalSingleScintillators(G4LogicalVolume *hcalenvelope);
  int CheckTiltAngle() const;
  int ConsistencyCheck() const;
  void SetTiltViaNcross();
  std::pair<int, int> GetLayerTowerId(G4VPhysicalVolume *volume) const;

 protected:
  int ConstructInnerHcal(G4LogicalVolume *sandwich);
  double x_at_y(Point_2 &p0, Point_2 &p1, double yin);
  PHG4InnerHcalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;
  G4AssemblyVolume *m_ScintiMotherAssembly = nullptr;
  double m_InnerRadius = NAN;
  double m_OuterRadius = NAN;
  double m_SizeZ = NAN;
  double m_ScintiTileX = NAN;
  double m_ScintiTileXLower = NAN;
  double m_ScintiTileXUpper = NAN;
  double m_ScintiTileZ = NAN;
  double m_ScintiTileThickness = NAN;
  double m_ScintiInnerGap = NAN;
  double m_ScintiOuterGap = NAN;
  double m_ScintiOuterRadius = NAN;
  double m_TiltAngle = NAN;
  double m_EnvelopeInnerRadius = NAN;
  double m_EnvelopeOuterRadius = NAN;
  double m_EnvelopeZ = NAN;
  double m_VolumeEnvelope = NAN;
  double m_VolumeSteel = NAN;
  double m_VolumeScintillator = NAN;

  int m_NumScintiPlates = 0;
  int m_NumScintiTilesPos = 0;
  int m_NumScintiTilesNeg = 0;

  int m_Active = 0;
  int m_AbsorberActive = 0;

  int m_Layer = 0;
  std::string m_SuperDetector;
  std::set<G4VPhysicalVolume *> m_SteelAbsorberPhysVolSet;
  std::map<G4VPhysicalVolume *, std::pair<int, int>> m_ScintiTilePhysVolMap;
  std::vector<G4VSolid *> m_ScintiTilesVec;
  std::string m_ScintiLogicNamePrefix;
};

#endif  // G4DETECTORS_PHG4INNERHCALDETECTOR_H
