#ifndef G4MVTX_PHG4MVTXSUPPORT_H
#define G4MVTX_PHG4MVTXSUPPORT_H

#include "ServiceStructure.h"
#include "Cable.h"

#include <PHG4MvtxDisplayAction.h>

#include <g4main/PHG4Detector.h>
//#include <g4detectors/PHG4DetectorGroupSubsystem.h>

#include <Geant4/G4AssemblyVolume.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4RotationMatrix.hh>  // for G4RotationMatrix
#include <Geant4/G4ThreeVector.hh>     // for G4ThreeVector
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4VSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

class PHG4MvtxSupport
{
 public:
  PHG4MvtxSupport(PHG4MvtxDisplayAction* dispAct);

  virtual ~PHG4MvtxSupport(){};

  void ConstructMvtxSupport(G4LogicalVolume *&lv);

 private:
  PHG4MvtxDisplayAction *m_DisplayAction;

  std::vector<float> get_thickness(ServiceStructure *object);
  void TrackingServiceCone(ServiceStructure *object, G4AssemblyVolume &assemblyVolume);
  void TrackingServiceCylinder(ServiceStructure *object, G4AssemblyVolume &assemblyVolume);
  void CreateCable(Cable *object, G4AssemblyVolume &assemblyVolume);
  void CreateCableBundle(G4AssemblyVolume &assemblyVolume, std::string superName, 
                         bool enableSignal, bool enableCooling, bool enablePower,
                         float x1, float x2, float y1, float y2, float z1, float z2);//, float theta = 0);
G4AssemblyVolume *buildBarrelCable();
G4AssemblyVolume *buildL0Cable();
G4AssemblyVolume *buildL1Cable();
G4AssemblyVolume *buildL2Cable();
G4AssemblyVolume *connectBarrelToLayer();
};

#endif
