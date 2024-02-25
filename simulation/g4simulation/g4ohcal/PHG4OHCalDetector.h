// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4OHCAL_PHG4OHCALDETECTOR_H
#define G4OHCAL_PHG4OHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <limits>
#include <map>
#include <set>
#include <string>  // for string
#include <tuple>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4OHCalDisplayAction;
class PHG4OHCalFieldSetup;
class PHParameters;
class PHG4Subsystem;
class PHG4GDMLConfig;
class RawTowerGeomContainer;

class PHG4OHCalDetector : public PHG4Detector
{
 public:
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
  std::tuple<int, int, int> GetRowColumnId(G4VPhysicalVolume *volume) const;

 private:
  void AddGeometryNode();
  int ConstructOHCal(G4LogicalVolume *hcalenvelope);
  int map_towerid(const int tower_id);
  int map_layerid(const unsigned int isector, const int layer_id);
  std::tuple<int, int, int> ExtractLayerTowerId(const unsigned int isector, G4VPhysicalVolume *volume);
  PHG4OHCalDisplayAction *m_DisplayAction{nullptr};
  PHG4OHCalFieldSetup *m_FieldSetup{nullptr};
  PHParameters *m_Params{nullptr};
  G4AssemblyVolume *m_ScintiMotherAssembly{nullptr};
  G4AssemblyVolume *m_ChimScintiMotherAssembly{nullptr};
  //! registry for volumes that should not be exported
  PHG4GDMLConfig *gdml_config{nullptr};
  RawTowerGeomContainer *m_RawTowerGeom{nullptr};

  double m_InnerRadius{std::numeric_limits<double>::quiet_NaN()};
  double m_OuterRadius{std::numeric_limits<double>::quiet_NaN()};
  double m_SizeZ{std::numeric_limits<double>::quiet_NaN()};
  double m_VolumeEnvelope{std::numeric_limits<double>::quiet_NaN()};
  double m_VolumeSteel{0};
  double m_VolumeScintillator{0};

  int m_NumScintiPlates{-9999};

  int m_ActiveFlag{0};
  int m_AbsorberActiveFlag{0};

  int m_Layer{0};
  std::string m_SuperDetector;
  std::set<G4LogicalVolume *> m_SteelAbsorberLogVolSet;
  std::set<G4LogicalVolume *> m_ScintiTileLogVolSet;
  std::map<G4VPhysicalVolume *, std::tuple<int, int, int>> m_ScintiTilePhysVolMap;

  std::string m_GDMPath;
  std::string m_TowerGeomNodeName;
};

#endif  // G4OHCAL_PHG4OHCALDETECTOR_H
