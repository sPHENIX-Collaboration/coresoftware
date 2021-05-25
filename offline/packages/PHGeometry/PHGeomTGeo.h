// $Id: $

/*!
 * \file PHGeomTGeo.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHGEOMETRY_PHGEOMTGEO_H
#define PHGEOMETRY_PHGEOMTGEO_H

#include <phool/PHObject.h>

#include <iostream>

class TGeoManager;

/*!
 * \brief PHGeomTGeo provide run-time access to TGeoManger.
 * It is transient object and it shall NOT be saved to DST.
 * For DST storage, please use PHGeomIOTGeo
 * For operation of this class with DST node, please use PHGeomUtility
 */
class PHGeomTGeo : public PHObject
{
 public:
  PHGeomTGeo();
  ~PHGeomTGeo() override;

  /** identify Function from PHObject
      @param os Output Stream
   */
  void identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  void Reset() override;

  /// isValid returns non zero if object contains vailid data
  int isValid() const override;

  //! Assign TGeoManager object.
  //! Once assigned, the TGeoManager will be locked to avoid a second TGeoManager override gGeoManager and lead to an invalid PHGeomTGeo
  //! The pointer TGeoManager should be the current gGeoManager
  void
  SetGeometry(TGeoManager* g);

  TGeoManager*
  GetGeometry();

  //! Since ROOT force TGeoManager is a unique object via global pointer gGeoManager,
  //! this function checks whether _fGeom is still the current gGeoManager which avoids operates on an invalid pointer
  bool ConsistencyCheck() const;

 protected:
  //! store and stream the full geometry via DST objects
  TGeoManager* _fGeom;
};

#endif
