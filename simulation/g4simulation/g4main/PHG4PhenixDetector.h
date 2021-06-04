// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4MAIN_PHG4PHENIXDETECTOR_H
#define G4MAIN_PHG4PHENIXDETECTOR_H

#include <Geant4/G4VUserDetectorConstruction.hh>
#include <Geant4/G4Types.hh>                      // for G4double

#include <list>
#include <string>                                 // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHG4Detector;
class PHG4PhenixDisplayAction;
class PHG4Reco;

//! this is the main detector construction class, passed to geant to construct the entire phenix detector
class PHG4PhenixDetector : public G4VUserDetectorConstruction
{
 public:
  //! constructor
  PHG4PhenixDetector(PHG4Reco* subsys);

  //! destructor
  ~PHG4PhenixDetector() override;

  void Verbosity(const int verb) { m_Verbosity = verb; }
  int Verbosity() const { return m_Verbosity; }

  //! register a detector. This is called in PHG4Reco::Init based on which detectors are found on the tree
  void AddDetector(PHG4Detector* detector)
  {
    m_DetectorList.push_back(detector);
  }

  //! this is called by geant to actually construct all detectors
  G4VPhysicalVolume* Construct() override;

  G4double GetWorldSizeX() const { return WorldSizeX; }

  G4double GetWorldSizeY() const { return WorldSizeY; }
  G4double GetWorldSizeZ() const { return WorldSizeZ; }

  void SetWorldSizeX(const G4double sx) { WorldSizeX = sx; }
  void SetWorldSizeY(const G4double sy) { WorldSizeY = sy; }
  void SetWorldSizeZ(const G4double sz) { WorldSizeZ = sz; }

  void SetWorldShape(const std::string& s) { worldshape = s; }
  void SetWorldMaterial(const std::string& s) { worldmaterial = s; }
  G4VPhysicalVolume* GetPhysicalVolume(void) { return physiWorld; }

 private:
  PHG4PhenixDisplayAction* m_DisplayAction;

  int m_Verbosity;

  //! list of detectors to be constructed

  std::list<PHG4Detector*> m_DetectorList;

  G4LogicalVolume* logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;  //pointer to the physical World
  G4double WorldSizeX;
  G4double WorldSizeY;
  G4double WorldSizeZ;
  std::string worldshape;
  std::string worldmaterial;
};

#endif  // G4MAIN_PHG4PHENIXDETECTOR_H
