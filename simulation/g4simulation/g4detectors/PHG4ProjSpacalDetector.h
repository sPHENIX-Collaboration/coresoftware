// $$Id: PHG4ProjSpacalDetector.h,v 1.2 2015/02/10 15:39:26 pinkenbu Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.2 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */

#ifndef PHG4ProjSpacalDetector_h
#define PHG4ProjSpacalDetector_h

#include "PHG4SpacalDetector.h"
#include "PHG4CylinderGeom_Spacalv2.h"

#include <Geant4/globals.hh>
#include <Geant4/G4Region.hh>
#include <Geant4/G4Types.hh>

#include <map>
#include <set>

class G4Material;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4UserLimits;

class PHG4ProjSpacalDetector : public PHG4SpacalDetector
{

public:
  typedef PHG4CylinderGeom_Spacalv2 SpacalGeom_t;

  PHG4ProjSpacalDetector(PHCompositeNode* Node, const std::string& dnam,
      SpacalGeom_t * geom, const int layer = 0);

  // empty dtor, step limits are deleted in base class
  virtual
    ~PHG4ProjSpacalDetector(void) {}

  virtual
  void
  Construct(G4LogicalVolume* world);

  virtual std::pair<G4LogicalVolume *, G4Transform3D>
  Construct_AzimuthalSeg();

  //! a block along z axis built with G4Trd that is slightly tapered in x dimension
 virtual G4LogicalVolume*
  Construct_Block();


  virtual
  void
  Print(const std::string& what = "ALL") const;

  const SpacalGeom_t *
  get_geom() const
  {
    return _geom;
  }

  virtual
  PHG4CylinderGeom * clone_geom() const
  {
    return new SpacalGeom_t(*_geom);
  }

private:

  SpacalGeom_t * _geom;

};

#endif
