#ifndef G4MVTX_PHG4MVTXSUPPORT_H
#define G4MVTX_PHG4MVTXSUPPORT_H

#include "PHG4MvtxDefs.h"

#include <array>
#include <string>

class G4AssemblyVolume;
class G4LogicalVolume;
class PHG4MvtxCable;
class PHG4MvtxDetector;
class PHG4MvtxDisplayAction;

class PHG4MvtxSupport
{
 public:
  PHG4MvtxSupport(PHG4MvtxDetector *detector, PHG4MvtxDisplayAction *dispAct, bool overlapCheck);

  virtual ~PHG4MvtxSupport();

  void ConstructMvtxSupport(G4LogicalVolume *&lv);

 private:
  PHG4MvtxDetector *m_Detector{nullptr};
  PHG4MvtxDisplayAction *m_DisplayAction{nullptr};

  void CreateMvtxSupportMaterials();

  void CreateEndWheelsSideN(G4AssemblyVolume *&av);
  void CreateEndWheelsSideS(G4AssemblyVolume *&av);
  void CreateConeLayers(G4AssemblyVolume *&av);
  void CreateCYSS(G4AssemblyVolume *&av);
  void CreateServiceBarrel(G4AssemblyVolume *&av);

  void GetEndWheelSideN(const int lay, G4AssemblyVolume *&endWheel);
  void GetEndWheelSideS(const int lay, G4AssemblyVolume *&endWheel);
  void GetConeVolume(int lay, G4AssemblyVolume *&av);

  void CreateCable(PHG4MvtxCable *object, G4AssemblyVolume &assemblyVolume);
  void CreateCableBundle(G4AssemblyVolume &assemblyVolume, const std::string &superName,
                         bool enableSignal, bool enableCooling, bool enablePower,
                         double x1, double x2, double y1, double y2, double z1, double z2);

  G4AssemblyVolume *buildBarrelCable();
  G4AssemblyVolume *buildLayerCables(const int &lay);

  G4AssemblyVolume *m_avSupport{nullptr};
  G4AssemblyVolume *m_avBarrelCable{nullptr};
  std::array<G4AssemblyVolume *, PHG4MvtxDefs::kNLayers> m_avLayerCable{nullptr, nullptr, nullptr};

  bool m_overlapCheck{false};
};

#endif
