// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4OUTERHCALDETECTOR_H
#define G4DETECTORS_PHG4OUTERHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>  // for G4double

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/point_generators_2.h>
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
class PHG4OHCalDisplayAction;
class PHG4OHCalFieldSetup;
class PHParameters;
class PHG4Subsystem;

class PHG4OHCalDetector : public PHG4Detector
{
 public:
  typedef CGAL::Exact_circular_kernel_2 Circular_k;
  typedef CGAL::Point_2<Circular_k> Point_2;
  //! constructor
  PHG4OHCalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam);

  //! destructor
  ~PHG4OHCalDetector() override;

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInOHCal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  int ConsistencyCheck() const;
  std::pair<int, int> GetLayerTowerId(G4VPhysicalVolume *volume) const;

 protected:

  int ConstructOHCal(G4LogicalVolume *hcalenvelope);
  PHG4OHCalDisplayAction *m_DisplayAction;
  PHG4OHCalFieldSetup *m_FieldSetup;
  PHParameters *m_Params;
  G4AssemblyVolume *m_ScintiMotherAssembly;
  G4AssemblyVolume *m_ChimScintiMotherAssembly;
  G4VSolid *m_SteelCutoutForMagnetG4Solid;
  double m_InnerRadius;
  double m_OuterRadius;
  double m_SizeZ;
  double m_ScintiTileX;
  double m_ScintiTileXLower;
  double m_ScintiTileXUpper;
  double m_ScintiTileZ;
  double m_ScintiTileThickness;
  double m_ScintiGap;
  double m_ScintiInnerRadius;
  double m_ScintiOuterRadius;
  double m_TiltAngle;
  double m_EnvelopeInnerRadius;
  double m_EnvelopeOuterRadius;
  double m_EnvelopeZ;
  double m_VolumeEnvelope;
  double m_VolumeSteel;
  double m_VolumeScintillator;

  int m_NumScintiPlates;
  int m_NumScintiTiles;

  int m_ActiveFlag;
  int m_AbsorberActiveFlag;

  int m_Layer;
  std::string m_SuperDetector;
  std::string m_ScintiLogicNamePrefix;
  std::vector<G4VSolid *> m_ScintiTilesVec;
  std::set<G4VPhysicalVolume *> m_SteelAbsorberVec;
  std::map<G4VPhysicalVolume *, std::pair<int, int>> m_ScintiTilePhysVolMap;

  std::string m_GDMPath;
};

#endif
