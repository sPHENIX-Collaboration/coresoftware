// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4OUTERHCALDETECTOR_H
#define G4DETECTORS_PHG4OUTERHCALDETECTOR_H

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
class PHG4OuterHcalDisplayAction;
class PHG4OuterHcalSubsystem;
class PHParameters;

class PHG4OuterHcalDetector : public PHG4Detector
{
 public:
  typedef CGAL::Exact_circular_kernel_2 Circular_k;
  typedef CGAL::Point_2<Circular_k> Point_2;
  //! constructor
  PHG4OuterHcalDetector(PHG4OuterHcalSubsystem *subsys, PHCompositeNode *Node, PHParameters *params, const std::string &dnam = "HCALOUT");

  //! destructor
  virtual ~PHG4OuterHcalDetector();

  //! construct
  virtual void Construct(G4LogicalVolume *world);

  virtual void Print(const std::string &what = "ALL") const;

  //!@name volume accessors
  //@{
  int IsInOuterHcal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
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
  PHG4OuterHcalDisplayAction *m_DisplayAction;
  PHG4OuterHcalFieldSetup *m_FieldSetup;
  PHParameters *m_Params;
  G4AssemblyVolume *m_ScintiMotherAssembly;
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

  int m_NScintiPlates;
  int m_NScintiTiles;

  int m_ActiveFlag;
  int m_AbsorberActiveFlag;

  int m_Layer;
  std::string m_SuperDetector;
  std::string m_ScintiLogicNamePrefix;
  std::vector<G4VSolid *> m_ScintiTilesVec;
  std::set<G4VPhysicalVolume *> m_SteelAbsorberVec;
};

#endif
