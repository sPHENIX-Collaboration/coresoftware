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
class PHG4ZDCDisplayAction;
class PHG4Subsystem;
class PHG4GDMLConfig;
class PHParameters;

/**
 * \file ${file_name}
 * \brief Module to build forward sampling Hadron calorimeterr (endcap) in Geant4
 * \author Nils Feege <nils.feege@stonybrook.edu>
 */

class PHG4ZDCDetector : public PHG4Detector
{
 public:
  //! constructor
  PHG4ZDCDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~PHG4ZDCDetector() {}

  //! construct
  virtual void ConstructMe(G4LogicalVolume *world);

  //!@name volume accessors
  int IsInZDC(G4VPhysicalVolume *) const;
  
  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  int get_Layer() const { return m_Layer; }

  PHG4ZDCDisplayAction *GetDisplayAction() { return m_DisplayAction; }

 private:
  G4LogicalVolume *ConstructTower(int type);
  

  PHG4ZDCDisplayAction *m_DisplayAction;
  PHParameters *m_Params;
  //! registry for volumes that should not be exported, i.e. fibers
  PHG4GDMLConfig *m_GdmlConfig;

  /* ZDC geometry */
  double m_Angle;

  double m_TPlate;
  double m_HPlate;
  double m_WPlate;
  
  double m_TAbsorber;
  double m_HAbsorber;
  double m_WAbsorber;

  double m_DFiber;
  double m_HFiber;
  double m_WFiber;
  double m_GFiber;

  double m_Gap;

  double m_XRot;
  double m_YRot;
  double m_ZRot;

  double m_PlaceX;
  double m_PlaceY;
  double m_PlaceZ;

  int m_NMod;
  int m_NLay;

  int m_ActiveFlag;
  int m_AbsorberActiveFlag;
  int m_Layer;

  std::string m_SuperDetector;
  

  std::set<G4LogicalVolume *> m_AbsorberLogicalVolSet;
  std::set<G4LogicalVolume *> m_ScintiLogicalVolSet;

 protected:
 
  PHParameters *GetParams() const { return m_Params; }
  void AbsorberLogicalVolSetInsert(G4LogicalVolume *logvol)
  {
    m_AbsorberLogicalVolSet.insert(logvol);
  }
  void ScintiLogicalVolSetInsert(G4LogicalVolume *logvol)
  {
    m_ScintiLogicalVolSet.insert(logvol);
  }
 
 
};

#endif
