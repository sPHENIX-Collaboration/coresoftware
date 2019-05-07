// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4MVTX_PHG4MVTXDETECTOR_H
#define G4MVTX_PHG4MVTXDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4RotationMatrix.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/globals.hh>

#include <array>
#include <map>
#include <set>
#include <string>
#include <vector>

class G4AssemblyVolume;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class PHG4MVTXDisplayAction;
class PHG4MVTXSubsystem;
class PHParameters;
class PHParametersContainer;

class PHG4MVTXDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4MVTXDetector(PHG4MVTXSubsystem* subsys, PHCompositeNode* Node, const PHParametersContainer* _paramsContainer, const std::string& dnam = "MVTX");

  //! destructor
  virtual ~PHG4MVTXDetector() {}

  //! construct
  virtual void Construct(G4LogicalVolume* world);

  //!@name volume accessors
  //@{
  int IsInMVTX(G4VPhysicalVolume*, int& layer, int& stave) const;
  int IsSensor(G4VPhysicalVolume*) const;
  //@}

  int IsActive(int lyr) const { return m_IsLayerActive[lyr]; }
  int IsAbsorberActive(int lyr) const { return m_IsLayerAbsorberActive[lyr]; }
  int IsBlackHole(int lyr) const { return m_IsBlackHole[lyr]; }
  void SuperDetector(const std::string& name) { superdetector = name; }
  const std::string SuperDetector() const { return superdetector; }
  void Detector(const std::string& name) { detector_type = name; }
  const std::string Detector() const { return detector_type; }

  int get_layer(int stv_index) const;
  int get_stave(int stv_index) const;

 private:
  void AddGeometryNode();
  int ConstructMVTX(G4LogicalVolume* sandwich);
  int ConstructMVTX_Layer(int layer, G4AssemblyVolume* stave, G4LogicalVolume*& trackerenvelope);
  void SetDisplayProperty(G4AssemblyVolume* av);
  void SetDisplayProperty(G4LogicalVolume* lv);
  void FillPVArray(G4AssemblyVolume* av);
  void FindSensor(G4LogicalVolume* lv);
  PHG4MVTXDisplayAction* m_DisplayAction;
  const PHParametersContainer* m_ParamsContainer;
  static constexpr int n_Layers = 3;

  // map of sensor physical volume pointers
  std::set<G4VPhysicalVolume*> m_SensorPV;
  std::map<G4VPhysicalVolume*, std::tuple<int, int>> m_StavePV;

  // setup parameters
  std::array<int, n_Layers> m_IsLayerActive;
  std::array<int, n_Layers> m_IsLayerAbsorberActive;
  std::array<int, n_Layers> m_IsBlackHole;
  std::array<int, n_Layers> m_N_staves;
  std::array<G4double, n_Layers> m_nominal_radius;
  std::array<G4double, n_Layers> m_nominal_phitilt;
  // sensor parameters
  double pixel_x;
  double pixel_z;
  double pixel_thickness;

  // calculated quantities
  G4double get_phistep(int lay) const { return 2.0 * M_PI / (double) m_N_staves[lay]; }

  std::string detector_type;
  std::string superdetector;
  std::string stave_geometry_file;
};

#endif
