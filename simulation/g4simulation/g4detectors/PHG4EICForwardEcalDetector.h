// Tell emacs that this is a C++ source
// -*- C++ -*-.
#ifndef G4DETECTORS_PHG4EICFORWARDECALDETECTOR_H
#define G4DETECTORS_PHG4EICFORWARDECALDETECTOR_H

#include <PHG4ForwardEcalDetector.h>

#include <Geant4/G4String.hh>
#include <Geant4/G4Types.hh>

#include <map>
#include <string>

class G4LogicalVolume;
class PHCompositeNode;
class PHG4ForwardEcalSubsystem;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4EICForwardEcalDetector : public PHG4ForwardEcalDetector
{
 public:
  //! constructor
  PHG4EICForwardEcalDetector(PHG4ForwardEcalSubsystem* subsys, PHCompositeNode* Node, const std::string& dnam = "BLOCK");

  //! destructor
  virtual ~PHG4EICForwardEcalDetector();

  //! construct
  virtual void Construct(G4LogicalVolume* world);

  void SetTowerDimensions(G4double dx, G4double dy, G4double dz)
  {
    _tower_dx = dx;
    _tower_dy = dy;
    _tower_dz = dz;
  }

  void SetMaterialScintillator(G4String material) { _materialScintillator = material; }
  void SetMaterialAbsorber(G4String material) { _materialAbsorber = material; }

 private:
  G4LogicalVolume* ConstructTower();
  int PlaceTower(G4LogicalVolume* envelope, G4LogicalVolume* tower);
  int ParseParametersFromTable();

  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
  };

  std::map<std::string, towerposition> _map_tower;

  /* ECAL tower geometry */
  G4double _tower_dx;
  G4double _tower_dy;
  G4double _tower_dz;

  G4String _materialScintillator;
  G4String _materialAbsorber;
};

#endif
