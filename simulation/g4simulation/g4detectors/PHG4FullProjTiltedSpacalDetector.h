// Tell emacs that this is a C++ source
//  -*- C++ -*-.

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */

#ifndef G4DETECTORS_PHG4FULLPROJTILTEDSPACALDETECTOR_H
#define G4DETECTORS_PHG4FULLPROJTILTEDSPACALDETECTOR_H

#include "PHG4CylinderGeom_Spacalv3.h"
#include "PHG4SpacalDetector.h"

#include <Geant4/G4Transform3D.hh>

#include <cassert>
#include <string>   // for string
#include <utility>  // for pair

class G4LogicalVolume;
class PHCompositeNode;
class PHG4CylinderGeom;
class PHG4Subsystem;
class PHParameters;

//! Fully projective SPACAL built from 2D tapered modules and allow azimuthal tilts
class PHG4FullProjTiltedSpacalDetector : public PHG4SpacalDetector
{
 public:
  typedef PHG4CylinderGeom_Spacalv3 SpacalGeom_t;

  PHG4FullProjTiltedSpacalDetector(PHG4Subsystem* subsys, PHCompositeNode* Node, const std::string& dnam,
                                   PHParameters* parameters, const int layer = 0);

  // empty dtor, step limits are deleted in base class
  ~PHG4FullProjTiltedSpacalDetector(void) override {}

  void
  ConstructMe(G4LogicalVolume* world) override;

  std::pair<G4LogicalVolume*, G4Transform3D>
  Construct_AzimuthalSeg() override;

  //! a block along z axis built with G4Trd that is slightly tapered in x dimension
  virtual G4LogicalVolume*
  Construct_Tower(const SpacalGeom_t::geom_tower& tower);
  //! a block for the light guide along z axis that fit to the tower
  virtual G4LogicalVolume*
  Construct_LightGuide(const SpacalGeom_t::geom_tower& tower, const int index_x, const int index_y);

  //! a block along z axis built with G4Trd that is slightly tapered in x dimension
  virtual int
  Construct_Fibers(const SpacalGeom_t::geom_tower& tower, G4LogicalVolume* LV_tower);

  //! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
  virtual int
  Construct_Fibers_SameLengthFiberPerTower(const SpacalGeom_t::geom_tower& tower, G4LogicalVolume* LV_tower);

  void
  Print(const std::string& what = "ALL") const override;

  PHG4CylinderGeom* clone_geom() const override
  {
    return new SpacalGeom_t(*get_geom_v3());
  }

 private:
  //  SpacalGeom_t* _geom;
  //! get the v3 cast of the geometry object
  SpacalGeom_t*
  get_geom_v3()
  {
    SpacalGeom_t* v3_geom = dynamic_cast<SpacalGeom_t*>(_geom);
    assert(v3_geom);
    return v3_geom;
  }

  const SpacalGeom_t*
  get_geom_v3() const
  {
    SpacalGeom_t* v3_geom = dynamic_cast<SpacalGeom_t*>(_geom);
    assert(v3_geom);
    return v3_geom;
  }
};

#endif
