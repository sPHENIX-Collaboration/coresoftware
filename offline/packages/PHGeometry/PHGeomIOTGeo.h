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
  PHGeomIOTGeo(const PHGeomIOTGeo& geom);
  virtual
  ~PHGeomIOTGeo();

  /// Virtual copy constructor.
  virtual PHObject*
  clone() const;

  /** identify Function from PHObject
   @param os Output Stream
   */
  virtual void
  identify(std::ostream& os = std::cout) const;

  /// Clear Event
  virtual void
  Reset();

  /// isValid returns non zero if object contains vailid data
  virtual int
  isValid() const;

  //! PHGeomIOTGeo do NOT own this TGeoVolume * g. Internally, it will use g to make a copy which PHGeomIOTGeo fully owns
  void
  SetGeometry(const TGeoVolume * g);

  //! Make a copy of TGeoVolume. The caller is responsible for deleting the returned TGeoVolume
  TGeoVolume *
  GetGeometryCopy() const;

  //! Get a constatnt copy of TGeoVolume. The caller is responsible for deleting the returned TGeoVolume
  const TGeoVolume *
  GetGeometry() const
  {
    return _fGeom;
  }

protected:

  //! store and stream the full geometry via DST objects
  TGeoVolume * _fGeom;

ClassDef(PHGeomIOTGeo,1)
};

#endif /* PHGeomIOTGeo_H_ */
