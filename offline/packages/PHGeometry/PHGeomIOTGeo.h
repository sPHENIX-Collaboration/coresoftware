// $Id: $

/*!
 * \file PHGeomIOTGeo.h
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef PHGEOMETRY_PHGEOMIOTGEO_H
#define PHGEOMETRY_PHGEOMIOTGEO_H

#include <phool/PHObject.h>

#include <iostream>
#include <vector>

class TGeoVolume;
class TGeoManager;

/*!
 * \brief PHGeomIOTGeo store geometry information to DST files in the format of binary streamed TGeoVolume. It completely owns the geometry object
 * For run-time use of TGeoManager, please use PHGeomTGeo
 * For operation of this class with DST node, please use PHGeomUtility
 */
class PHGeomIOTGeo : public PHObject
{
 public:
  PHGeomIOTGeo();
  ~PHGeomIOTGeo() override;

  PHObject* CloneMe() const override {return new PHGeomIOTGeo(*this);}

  /** identify Function from PHObject
   @param os Output Stream
   */
  void
  identify(std::ostream& os = std::cout) const override;

  /// Clear Event
  void
  Reset() override;

  /// isValid returns non zero if object contains vailid data
  int
  isValid() const override;

  //! PHGeomIOTGeo do NOT own this TGeoVolume * g. Internally, it will use g to make a copy which PHGeomIOTGeo fully owns
  void
  SetGeometry(const TGeoVolume* g);

  //! Construct TGeoManager. The result TGeoManager is not yet closed and open for further editing
  TGeoManager*
  ConstructTGeoManager();

  //! Make a copy of TGeoVolume.
  //! The caller is responsible for deleting the returned TGeoVolume
  //! The caller is also responsible for constructing a valid TGeoManager before calling this function
  TGeoVolume*
  GetGeometryCopy();

  std::vector<char>&
  GetData()
  {
    return Data;
  }

  const std::vector<char>&
  GetData() const
  {
    return Data;
  }

 protected:
  //! store the streamed geometry and its streamer via a binary stream
  std::vector<char> Data;

  ClassDefOverride(PHGeomIOTGeo, 3)
};

#endif
