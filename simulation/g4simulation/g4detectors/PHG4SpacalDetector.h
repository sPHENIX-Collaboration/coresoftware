// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2014/08/12 03:49:12 $$
 */

#ifndef G4DETECTORS_PHG4SPACALDETECTOR_H
#define G4DETECTORS_PHG4SPACALDETECTOR_H

#include "PHG4CylinderGeom_Spacalv1.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Transform3D.hh>
#include <Geant4/G4Types.hh>  // for G4double

#include <map>
#include <string>   // for string
#include <utility>  // for pair

class G4LogicalVolume;
class G4UserLimits;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4CylinderGeom;
class PHG4GDMLConfig;
class PHG4SpacalDisplayAction;
class PHParameters;
class PHG4Subsystem;

class PHG4SpacalDetector : public PHG4Detector
{
 public:
  typedef PHG4CylinderGeom_Spacalv1 SpacalGeom_t;

  PHG4SpacalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam,
                     PHParameters* parameters, const int layer = 0, bool init_geom = true);

  ~PHG4SpacalDetector(void) override;

  void
  ConstructMe(G4LogicalVolume* world) override;

  virtual std::pair<G4LogicalVolume*, G4Transform3D>
  Construct_AzimuthalSeg();

  virtual G4LogicalVolume*
  Construct_Fiber(const G4double length, const std::string& id);

  void
  SetActive(const int i = 1)
  {
    active = i;
  }

  void
  SetAbsorberActive(const int i = 1)
  {
    absorberactive = i;
  }

  void
  SetDetectorType(const std::string& typ)
  {
    detector_type = typ;
  }

  int IsInCylinderActive(const G4VPhysicalVolume*);

  void
  SuperDetector(const std::string& name)
  {
    superdetector = name;
  }

  const std::string
  SuperDetector() const
  {
    return superdetector;
  }

  int get_Layer() const
  {
    return layer;
  }

  void
  Print(const std::string& what = "ALL") const override;

  const SpacalGeom_t*
  get_geom() const
  {
    return _geom;
  }

  virtual PHG4CylinderGeom* clone_geom() const
  {
    return new SpacalGeom_t(*_geom);
  }

  enum
  {
    FIBER_CORE = 1,
    FIBER_CLADING = 0,
    ABSORBER = -1,
    SUPPORT = -2,
    INACTIVE = -100
  };

  PHG4SpacalDisplayAction* GetDisplayAction() { return m_DisplayAction; }

  void CosmicSetup(const int i) { m_CosmicSetupFlag = i; }
  int CosmicSetup() { return m_CosmicSetupFlag; }

 private:
  PHG4SpacalDisplayAction* m_DisplayAction = nullptr;

 protected:
  std::map<const G4VPhysicalVolume*, int> fiber_core_vol;

  //! map for G4VPhysicalVolume -> fiber ID
  std::map<const G4VPhysicalVolume*, int> fiber_vol;

  //! map for G4VPhysicalVolume -> Sector ID
  std::map<const G4VPhysicalVolume*, int> calo_vol;

  //! map for G4VPhysicalVolume -> towers ID
  std::map<const G4VPhysicalVolume*, int> block_vol;

  int active = 0;
  int absorberactive = 0;
  int layer = -9999;
  int m_CosmicSetupFlag = 0;
  std::string detector_type;
  std::string superdetector;

  //  G4UserLimits * step_limits;
  //  G4UserLimits * clading_step_limits;
  G4UserLimits* fiber_core_step_limits = nullptr;

  //! registry for volumes that should not be exported, i.e. fibers
  PHG4GDMLConfig* gdml_config = nullptr;
  //private:

  SpacalGeom_t* _geom = nullptr;
};

#endif
