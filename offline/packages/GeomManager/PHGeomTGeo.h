// $Id: $                                                                                             

/*!
 * \file PHGeomTGeo.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHGEOMTGEO_H_
#define PHGEOMTGEO_H_

#include <phool/PHObject.h>

class TGeoManager;

/*!
 * \brief PHGeomTGeo
 */
class PHGeomTGeo : public PHObject
{
public:
  PHGeomTGeo();
  virtual
  ~PHGeomTGeo();

  /** identify Function from PHObject
      @param os Output Stream
   */
  virtual void identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void Reset();

  /// isValid returns non zero if object contains vailid data
  virtual int isValid() const;

  //! WARNING : the new pointer do not belong to this class,
  //! the TGeoManager should be only deleted after PdbFvtxAlignment is deleted or set to another Geometry object
  void
  SetGeometry(TGeoManager * g);

  TGeoManager *
  GetGeometry();

protected:

  //! store and stream the full geometry via DST objects
  TGeoManager * _fGeom;

  ClassDef(PHGeomTGeo,1)
};

#endif /* PHGEOMTGEO_H_ */
