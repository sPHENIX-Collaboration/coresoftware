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
    int idx_j;
    int idx_k;
  };

  PHG4ForwardHcalDisplayAction *m_DisplayAction = nullptr;
  PHParameters *m_Params = nullptr;

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
