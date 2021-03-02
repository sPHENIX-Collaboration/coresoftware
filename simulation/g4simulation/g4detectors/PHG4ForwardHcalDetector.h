// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDHCALDETECTOR_H
#define G4DETECTORS_PHG4FORWARDHCALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>   // for G4double

#include <map>
#include <set>
#include <string>

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ForwardHcalDisplayAction;
class PHG4Subsystem;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4ForwardHcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4ForwardHcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4ForwardHcalDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInForwardHcal(G4VPhysicalVolume *) const;

  void SetXRot(G4double rot_in_x) { m_XRot = rot_in_x; }
  void SetYRot(G4double rot_in_y) { m_YRot = rot_in_y; }
  void SetZRot(G4double rot_in_z) { m_ZRot = rot_in_z; }

  void SetMaterialScintillator(const std::string &material) { m_MaterialScintillator = material; }
  void SetMaterialAbsorber(const std::string &material) { m_MaterialAbsorber = material; }

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  int get_Layer() const { return m_Layer; }

 private:
  G4LogicalVolume *ConstructTower();
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower);
  int ParseParametersFromTable();

  struct towerposition
  {
    G4double x;
    G4double y;
    G4double z;
  };

  PHG4ForwardHcalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;

  /* Calorimeter envelope geometry */
  G4double m_PlaceX;
  G4double m_PlaceY;
  G4double m_PlaceZ;

  G4double m_XRot;
  G4double m_YRot;
  G4double m_ZRot;

  G4double m_RMin1;
  G4double m_RMax1;
  G4double m_RMin2;
  G4double m_RMax2;

  G4double m_dZ;
  G4double m_SPhi;
  G4double m_DPhi;

  G4double m_WlsDw;
  G4double m_SupportDw;

  std::string m_MaterialScintillator;
  std::string m_MaterialAbsorber;

  int m_ActiveFlag = 1;
  int m_AbsorberActiveFlag = 0;
  int m_Layer = 0;

  std::string m_TowerLogicNamePrefix;
  std::string m_SuperDetector;

  std::map<std::string, G4double> m_GlobalParameterMap;
  std::map<std::string, towerposition> m_TowerPostionMap;

  std::set<G4LogicalVolume *> m_AbsorberLogicalVolSet;
  std::set<G4LogicalVolume *> m_ScintiLogicalVolSet;
};

#endif
