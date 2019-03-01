// $$Id: PHG4FullProjSpacalDetector.h,v 1.2 2015/02/10 15:39:26 pinkenbu Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */

#ifndef PHG4FullProjSpacalDetector_h
#define PHG4FullProjSpacalDetector_h

#include "PHG4CylinderGeom_Spacalv3.h"
#include "PHG4SpacalDetector.h"

#include <Geant4/G4Region.hh>
#include <Geant4/G4Types.hh>
#include <Geant4/globals.hh>

#include <cassert>
#include <map>
#include <set>

class G4Material;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;
class PHParameters;

//! Fully projective SPACAL built from 2D tapered modules.
//! This class is obsolete and for comparison study only. Use PHG4FullProjTiltedSpacalDetector instead.
//! It loads Chris Cullen 2D spacal design July 2015 by default.
class PHG4FullProjSpacalDetector : public PHG4SpacalDetector
{
 public:
  typedef PHG4CylinderGeom_Spacalv3 SpacalGeom_t;

  PHG4FullProjSpacalDetector(PHCompositeNode* Node, const std::string& dnam,
                             PHParameters* parameters, const int layer = 0);

  // empty dtor, step limits are deleted in base class
  virtual ~PHG4FullProjSpacalDetector(void) {}
  virtual void
  Construct(G4LogicalVolume* world);

  virtual std::pair<G4LogicalVolume*, G4Transform3D>
  Construct_AzimuthalSeg();

  //! a block along z axis built with G4Trd that is slightly tapered in x dimension
  virtual G4LogicalVolume*
  Construct_Tower(const SpacalGeom_t::geom_tower& tower);

  //! a block along z axis built with G4Trd that is slightly tapered in x dimension
  virtual int
  Construct_Fibers(const SpacalGeom_t::geom_tower& tower, G4LogicalVolume* LV_tower);

  //! Fully projective spacal with 2D tapered modules. To speed up construction, same-length fiber is used cross one tower
  virtual int
  Construct_Fibers_SameLengthFiberPerTower(const SpacalGeom_t::geom_tower& tower, G4LogicalVolume* LV_tower);

  virtual void
  Print(const std::string& what = "ALL") const;

  virtual PHG4CylinderGeom* clone_geom() const
  {
    return new SpacalGeom_t(*get_geom_v3());
  }

 private:
//  SpacalGeom_t* _geom;

  //! get the v3 cast of the geometry object
  SpacalGeom_t *
  get_geom_v3()
  {
    SpacalGeom_t *v3_geom =  dynamic_cast<SpacalGeom_t *> (_geom);
    assert(v3_geom);
    return v3_geom;
  }

  const SpacalGeom_t *
  get_geom_v3() const
  {
    SpacalGeom_t *v3_geom =  dynamic_cast<SpacalGeom_t *> (_geom);
    assert(v3_geom);
    return v3_geom;
  }


};

#endif
