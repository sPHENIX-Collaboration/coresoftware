#ifndef G4MVTX_PHG4MVTXSUPPORT_H
#define G4MVTX_PHG4MVTXSUPPORT_H

#include "PHG4MvtxDefs.h"

#include <Geant4/G4AssemblyVolume.hh>

#include <string>
#include <vector>
#include <array>

class G4LogicalVolume;
class PHG4MvtxCable;
class PHG4MvtxDisplayAction;
class PHG4MvtxServiceStructure;

class PHG4MvtxSupport
{
 public:
  PHG4MvtxSupport( PHG4MvtxDisplayAction *dispAct, bool overlapCheck );

  ~PHG4MvtxSupport();

  void ConstructMvtxSupport( G4LogicalVolume *&lv );

 private:
  PHG4MvtxDisplayAction *m_DisplayAction;

  void CreateMvtxSupportMaterials();

  void CreateEndWheelsSideN( G4AssemblyVolume *&av );
  void CreateEndWheelsSideS( G4AssemblyVolume *&av );
  void CreateConeLayers( G4AssemblyVolume *&av );
  void CreateCYSS( G4AssemblyVolume *&av );
  void CreateServiceBarrel( G4AssemblyVolume *&av );

  void GetEndWheelSideN( const int lay, G4AssemblyVolume *&endWheel );
  void GetEndWheelSideS( const int lay, G4AssemblyVolume *&endWheel );
  void GetConeVolume( int lay, G4AssemblyVolume *& av );

  void CreateCable( PHG4MvtxCable *object, G4AssemblyVolume &assemblyVolume );
  void CreateCableBundle( G4AssemblyVolume &assemblyVolume, const std::string &superName,
                          bool enableSignal, bool enableCooling, bool enablePower,
                          float x1, float x2, float y1, float y2, float z1, float z2);

  G4AssemblyVolume *buildBarrelCable();
  G4AssemblyVolume *buildLayerCables( const int &lay);

  G4AssemblyVolume *m_avSupport;
  G4AssemblyVolume *m_avBarrelCable;
  std::array<G4AssemblyVolume*, PHG4MvtxDefs::kNLayers> m_avLayerCable;

  bool m_overlapCheck = false;
};

#endif
