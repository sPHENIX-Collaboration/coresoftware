// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4PROJCRYSTALCALORIMETERDETECTOR_H
#define G4DETECTORS_PHG4PROJCRYSTALCALORIMETERDETECTOR_H

#include "PHG4CrystalCalorimeterDetector.h"

#include "PHG4CrystalCalorimeterDefs.h"

#include <Geant4/G4Types.hh>  // for G4double, G4int

#include <string>  // for string

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4Subsystem;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build crystal calorimeter with projective crystal geometry in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4ProjCrystalCalorimeterDetector : public PHG4CrystalCalorimeterDetector
{
 public:
  //! constructor
  PHG4ProjCrystalCalorimeterDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, PHParameters* parameters, const std::string& dnam);

  //! destructor
  virtual ~PHG4ProjCrystalCalorimeterDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume* world) override;

  //!@name volume accessors
  virtual int IsInCrystalCalorimeter(G4VPhysicalVolume*) const override;

  // ----- accessing member variables: ------------

  void SetCrystalSize(G4double dx_front, G4double dy_front, G4double dx_back, G4double dy_back, G4double dz)
  {
    _dx_front = dx_front;
    _dy_front = dy_front;
    _dx_back = dx_back;
    _dy_back = dy_back;
    _dz_crystal = dz;
  }

  void GetCrystalSize(G4double& dx_front, G4double& dy_front, G4double& dx_back, G4double& dy_back, G4double& dz)
  {
    dx_front = _dx_front;
    dy_front = _dy_front;
    dx_back = _dx_back;
    dy_back = _dy_back;
    dz = _dz_crystal;
  }

  void GetCarbonFiberAdjustments(G4double& adjust_width, G4double& adjust_length);

  void GetCarbonFiberSpacing(G4double& CF_width, G4double& Air_CF, G4double& Air_Cry);

protected:
  int GetCaloType() const override {return PHG4CrystalCalorimeterDefs::CaloType::projective;}

 private:
  int ConstructProjectiveCrystals(G4LogicalVolume* envelope);
  int Fill4x4Unit(G4LogicalVolume* crystal_logic);
  int FillSpecialUnit(G4LogicalVolume* crystal_logic, G4int ident);

  /* crystal geometry */
  G4double _dx_front;
  G4double _dy_front;
  G4double _dx_back;
  G4double _dy_back;
  G4double _dz_crystal;

  std::string _crystallogicnameprefix;

  bool _overlapcheck_local;
  std::set<G4VPhysicalVolume*> m_ActiveVolumeSet;
  std::set<G4VPhysicalVolume*> m_PassiveVolumeSet;
  // since getting parameters is a map search we do not want to
  // do this in every step, the parameters used are cached
  // in the following variables
  int m_IsActive;
  int m_AbsorberActive;
};

#endif
