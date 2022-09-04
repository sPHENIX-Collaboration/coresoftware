#ifndef G4MVTX_PHG4MVTXSUPPORT_H
#define G4MVTX_PHG4MVTXSUPPORT_H

#include <string>
#include <vector>

class G4AssemblyVolume;
class G4LogicalVolume;
class PHG4MvtxCable;
class PHG4MvtxDisplayAction;
class PHG4MvtxServiceStructure;

class PHG4MvtxSupport
{
 public:
  PHG4MvtxSupport(PHG4MvtxDisplayAction *dispAct, bool overlapCheck);

  ~PHG4MvtxSupport();

  void ConstructMvtxSupport(G4LogicalVolume *&lv);

 private:
  PHG4MvtxDisplayAction *m_DisplayAction;

  std::vector<float> get_thickness(PHG4MvtxServiceStructure *object);
  void TrackingServiceCone(PHG4MvtxServiceStructure *object, G4AssemblyVolume &assemblyVolume);
  void TrackingServiceCylinder(PHG4MvtxServiceStructure *object, G4AssemblyVolume &assemblyVolume);
  void CreateCable(PHG4MvtxCable *object, G4AssemblyVolume &assemblyVolume);
  void CreateCableBundle(G4AssemblyVolume &assemblyVolume, const std::string &superName,
                         bool enableSignal, bool enableCooling, bool enablePower,
                         float x1, float x2, float y1, float y2, float z1, float z2);  //, float theta = 0);

  G4AssemblyVolume *buildBarrelCable();
  G4AssemblyVolume *buildL0Cable();
  G4AssemblyVolume *buildL1Cable();
  G4AssemblyVolume *buildL2Cable();

  G4AssemblyVolume *m_avSupport;
  G4AssemblyVolume *m_avBarrelCable;
  G4AssemblyVolume *m_avL0Cable;
  G4AssemblyVolume *m_avL1Cable;
  G4AssemblyVolume *m_avL2Cable;
  //std::vector<G4AssemblyVolume*> m_endwheelCable;

  bool m_overlapCheck = false;
};

#endif
