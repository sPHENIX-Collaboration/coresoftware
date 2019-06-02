// Tell emacs that this is a C++ source
//  -*- C++ -*-.

#ifndef G4DETECTORS_PHG4SPACALPROTOTYPEDETECTOR_H
#define G4DETECTORS_PHG4SPACALPROTOTYPEDETECTOR_H

#include "PHG4CylinderGeom_Spacalv3.h"

#include <g4main/PHG4Detector.h>

#include <Geant4/G4Types.hh>
#include <Geant4/G4Transform3D.hh>

#include <map>
#include <string>                       // for string
#include <utility>

class G4LogicalVolume;
class G4UserLimits;
class G4VPhysicalVolume;
class G4VSolid;
class PHCompositeNode;
class PHParameters;

class PHG4SpacalPrototypeDetector : public PHG4Detector
{

public:
  typedef PHG4CylinderGeom_Spacalv3 SpacalGeom_t;

  PHG4SpacalPrototypeDetector(PHCompositeNode* Node, PHParameters *parameters, const std::string& dnam);

  virtual
  ~PHG4SpacalPrototypeDetector(void);

  virtual
  void
  Construct(G4LogicalVolume* world);

  virtual std::pair<G4LogicalVolume *, G4Transform3D>
  Construct_AzimuthalSeg();

  //! a block along z axis built with G4Trd that is slightly tapered in x dimension
  virtual G4LogicalVolume*
  Construct_Tower(const SpacalGeom_t::geom_tower & tower);
  //! a block for the light guide along z axis that fit to the tower
  virtual G4LogicalVolume*
  Construct_LightGuide(const SpacalGeom_t::geom_tower & tower, const int index_x, const int index_y);

  //! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
  virtual int
  Construct_Fibers_SameLengthFiberPerTower(
      const SpacalGeom_t::geom_tower & tower, G4LogicalVolume* LV_tower);

  virtual
  G4LogicalVolume *
  Construct_Fiber(const G4double length, const std::string & id);

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

  int
  IsInCylinderActive(const G4VPhysicalVolume*);

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

  virtual void
  Print(const std::string& what = "ALL") const;

  const SpacalGeom_t *
  get_geom() const
  {
    return _geom;
  }

  enum
  {
    FIBER_CORE = 1,
    FIBER_CLADING = 0,
    ABSORBER = -1,
    SUPPORT = -2,
    INACTIVE = -100
  };

protected:

  PHParameters *construction_params;

  G4VSolid* cylinder_solid;
  G4LogicalVolume* cylinder_logic;
  G4VPhysicalVolume* cylinder_physi;
  std::map<const G4VPhysicalVolume*, int> fiber_core_vol;

  //! map for G4VPhysicalVolume -> fiber ID
  std::map<const G4VPhysicalVolume*, int> fiber_vol;

  //! map for G4VPhysicalVolume -> Sector ID
  std::map<const G4VPhysicalVolume*, int> calo_vol;

  //! map for G4VPhysicalVolume -> towers ID
  std::map<const G4VPhysicalVolume*, int> block_vol;

  int active;
  int absorberactive;
  std::string detector_type;
  std::string superdetector;

  G4UserLimits * step_limits;
  G4UserLimits * clading_step_limits;
  G4UserLimits * fiber_core_step_limits;

private:
  SpacalGeom_t * _geom;

};

#endif
