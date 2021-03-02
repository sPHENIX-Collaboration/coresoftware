// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef G4DETECTORS_PHG4FORWARDECALDETECTOR_H
#define G4DETECTORS_PHG4FORWARDECALDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <cassert>
#include <map>
#include <set>
#include <string>
#include <utility>                // for pair, make_pair

class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4ForwardEcalDisplayAction;
class PHG4Subsystem;
class PHG4GDMLConfig;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4ForwardEcalDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4ForwardEcalDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4ForwardEcalDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInForwardEcal(G4VPhysicalVolume *) const;

  void SetTowerDimensions(double dx, double dy, double dz, int type);

  void SetPlace(double place_in_x, double place_in_y, double place_in_z)
  {
    m_PlaceX = place_in_x;
    m_PlaceY = place_in_y;
    m_PlaceZ = place_in_z;
  }
  void SetPlaceX(double place_x) { m_PlaceX = place_x; }
  void SetPlaceY(double place_y) { m_PlaceY = place_y; }
  void SetPlaceZ(double place_z) { m_PlaceZ = place_z; }

  double GetPlaceX() const { return m_PlaceX; }
  double GetPlaceY() const { return m_PlaceY; }
  double GetPlaceZ() const { return m_PlaceZ; }

  void SetXRot(double rot_in_x) { m_XRot = rot_in_x; }
  void SetYRot(double rot_in_y) { m_YRot = rot_in_y; }
  void SetZRot(double rot_in_z) { m_ZRot = rot_in_z; }

  double GetXRot() const { return m_XRot; }
  double GetYRot() const { return m_YRot; }
  double GetZRot() const { return m_ZRot; }

  double GetRMin(int i) const
  {
    assert(i >= 0 && i <= 1);
    return m_RMin[i];
  }
  double GetRMax(int i) const
  {
    assert(i >= 0 && i <= 1);
    return m_RMax[i];
  }

  void SetRMin(int i, double val)
  {
    assert(i >= 0 && i <= 1);
    m_RMin[i] = val;
  }
  void SetRMax(int i, double val)
  {
    assert(i >= 0 && i <= 1);
    m_RMax[i] = val;
  }

  double GetdZ() const { return m_dZ; }

  void SetdZ(double val) { m_dZ = val; }

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  int get_Layer() const { return m_Layer; }

  PHG4ForwardEcalDisplayAction *GetDisplayAction() { return m_DisplayAction; }

 private:
  G4LogicalVolume *ConstructTower(int type);
  G4LogicalVolume *ConstructTowerType2();
  G4LogicalVolume *ConstructTowerType3_4_5_6(int type);
  int PlaceTower(G4LogicalVolume *envelope, G4LogicalVolume *tower[6]);
  int ParseParametersFromTable();

  struct towerposition
  {
    double x;
    double y;
    double z;
    int type;
    int idx_j;
    int idx_k;
  };

  PHG4ForwardEcalDisplayAction *m_DisplayAction;
  PHParameters *m_Params;
  //! registry for volumes that should not be exported, i.e. fibers
  PHG4GDMLConfig *m_GdmlConfig;

  /* ECAL tower geometry */
  double m_TowerDx[7];
  double m_TowerDy[7];
  double m_TowerDz[7];

  double m_XRot;
  double m_YRot;
  double m_ZRot;

  double m_PlaceX;
  double m_PlaceY;
  double m_PlaceZ;

  double m_RMin[2];
  double m_RMax[2];

  double m_dZ;

  int m_ActiveFlag;
  int m_AbsorberActiveFlag;
  int m_Layer;

  std::string m_SuperDetector;
  std::string m_TowerLogicNamePrefix;

  std::map<std::string, towerposition> m_TowerPositionMap;
  std::map<std::string, double> m_GlobalParameterMap;

  std::set<G4LogicalVolume *> m_AbsorberLogicalVolSet;
  std::set<G4LogicalVolume *> m_ScintiLogicalVolSet;

 protected:
  const std::string TowerLogicNamePrefix() const { return m_TowerLogicNamePrefix; }
  PHParameters *GetParams() const { return m_Params; }
  void AbsorberLogicalVolSetInsert(G4LogicalVolume *logvol)
  {
    m_AbsorberLogicalVolSet.insert(logvol);
  }
  void ScintiLogicalVolSetInsert(G4LogicalVolume *logvol)
  {
    m_ScintiLogicalVolSet.insert(logvol);
  }
  std::map<std::string, double>::const_iterator FindIter(const std::string &name) { return m_GlobalParameterMap.find(name); }
  std::map<std::string, double>::const_iterator EndIter() { return m_GlobalParameterMap.end(); }
  void InsertParam(const std::string &parname, double parval) { m_GlobalParameterMap.insert(make_pair(parname, parval)); }
};

#endif
