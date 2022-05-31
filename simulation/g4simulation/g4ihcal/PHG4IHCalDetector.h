// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4INNERHCALDETECTOR_H
#define G4DETECTORS_PHG4INNERHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cmath>
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
class PHG4IHCalDisplayAction;
class PHParameters;
class PHG4Subsystem;

class PHG4IHCalDetector : public PHG4Detector
{
 public:

  //! constructor
  PHG4IHCalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  ~PHG4IHCalDetector() override;

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInIHCal(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }
  int get_Layer() const { return m_Layer; }
  G4AssemblyVolume *ConstructHcalScintillatorAssembly(G4LogicalVolume *hcalenvelope);
  void ConstructHcalSingleScintillators(G4LogicalVolume *hcalenvelope);
  int ConsistencyCheck() const;
  std::pair<int, int> GetLayerTowerId(G4VPhysicalVolume *volume) const;

 protected:

  int ConstructAbsorber(G4AssemblyVolume *avol, G4LogicalVolume *hcalenvelope);
  int ConstructScinTiles(G4AssemblyVolume *avol, G4LogicalVolume *hcalenvelope);

  int ConstructIHCal(G4LogicalVolume *sandwich);
  PHG4IHCalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;
  G4AssemblyVolume *m_ScintiMotherAssembly = nullptr;
  double m_InnerRadius = NAN;
  double m_OuterRadius = NAN;
  double m_SizeZ = NAN;
  double m_ScintiTileZ = NAN;
/*
  double m_ScintiTileXLower = NAN;
  double m_ScintiTileXUpper = NAN;
  double m_ScintiTileThickness = NAN;
  double m_ScintiInnerGap = NAN;
  double m_ScintiOuterGap = NAN;
  double m_ScintiOuterRadius = NAN;
  double m_TiltAngle = NAN;
*/
  double m_EnvelopeInnerRadius = NAN;
  double m_EnvelopeOuterRadius = NAN;
  double m_EnvelopeZ = NAN;
  double m_VolumeEnvelope = NAN;
  double m_VolumeSteel = NAN;
  double m_VolumeScintillator = NAN;

  int m_NumScintiPlates = -9999;
  int m_NumScintiTilesPos = -9999;
  int m_NumScintiTilesNeg = -9999;

  int m_Active = 0;
  int m_AbsorberActive = 0;

  int m_Layer = 0;
  std::string m_SuperDetector;
  std::set<G4VPhysicalVolume *> m_SteelAbsorberPhysVolSet;
  std::map<G4VPhysicalVolume *, std::pair<int, int>> m_ScintiTilePhysVolMap;
  std::vector<G4VSolid *> m_ScintiTilesVec;
  std::string m_ScintiLogicNamePrefix = "HcalInnerScinti";

  std::string m_GDMPath;
};

#endif  // G4DETECTORS_PHG4INNERHCALDETECTOR_H
