// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4IHCAL_PHG4IHCALDETECTOR_H
#define G4IHCAL_PHG4IHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cmath>
#include <map>
#include <set>
#include <string>  // for string
#include <tuple>
#include <utility>  // for pair

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4IHCalDisplayAction;
class PHParameters;
class PHG4Subsystem;
class PHG4GDMLConfig;

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
  std::tuple<int, int, int> GetLayerTowerId(G4VPhysicalVolume *volume) const;
  int GetSectorId(G4VPhysicalVolume *volume) const;

 private:
  int map_towerid(const int tower_id);
  int map_layerid(const int layer_id);
  int ConstructIHCal(G4LogicalVolume *sandwich);
  std::tuple<int, int, int> ExtractLayerTowerId(const unsigned int isector, G4VPhysicalVolume *volume);
  PHG4IHCalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;
  G4AssemblyVolume *m_ScintiMotherAssembly = nullptr;
  double m_InnerRadius = NAN;
  double m_OuterRadius = NAN;
  double m_SizeZ = NAN;
  double m_VolumeEnvelope = NAN;
  double m_VolumeSteel = 0.;
  double m_VolumeScintillator = 0;

  int m_NumScintiPlates = -9999;

  int m_Active = 0;
  int m_AbsorberActive = 0;

  int m_Layer = 0;

  std::string m_SuperDetector;
  std::set<G4LogicalVolume *> m_SteelAbsorberLogVolSet;
  std::set<G4LogicalVolume *> m_ScintiTileLogVolSet;
  std::map<G4VPhysicalVolume *, std::tuple<int, int, int>> m_ScintiTilePhysVolMap;
  std::map<G4VPhysicalVolume *, int> m_AbsorberPhysVolMap;

  //! registry for volumes that should not be exported
  PHG4GDMLConfig *gdml_config = nullptr;

  std::string m_GDMPath;
};

#endif  // G4IHCAL_PHG4IHCALDETECTOR_H
