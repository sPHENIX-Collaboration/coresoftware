// $Id: $                                                                                             

/*!
 * \file PHGeomIOTGeo.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHGeomIOTGeo_H_
#define PHGeomIOTGeo_H_

#include <phool/PHObject.h>

class TGeoVolume;

/*!
 * \brief PHGeomIOTGeo store geometry information to DST files in the format of TGeoVolume. It completely owns the geometry object
 */
class PHGeomIOTGeo : public PHObject
{
public:
  PHGeomIOTGeo();
  virtual
  ~PHGeomIOTGeo();

  /** identify Function from PHObject
      @param os Output Stream
   */
  virtual void identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void Reset();

  /// isValid returns non zero if object contains vailid data
  virtual int isValid() const;

  //! WARNING : the new pointer do not belong to this class,
  //! the TGeoVolume should be only deleted after PdbFvtxAlignment is deleted or set to another Geometry object
  void
  SetGeometry(TGeoVolume * g);

  TGeoVolume *
  GetGeometry();

protected:

  //! store and stream the full geometry via DST objects
  TGeoVolume * _fGeom;

  ClassDef(PHGeomIOTGeo,1)
};

#endif /* PHGeomIOTGeo_H_ */
