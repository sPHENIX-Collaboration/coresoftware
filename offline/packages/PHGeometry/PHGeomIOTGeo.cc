// $Id: $                                                                                             

/*!
 * \file PHGeomIOTGeo.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHGeomIOTGeo.h"

#include "TGeoVolume.h"
#include <iostream>

using namespace std;

ClassImp(PHGeomIOTGeo);

PHGeomIOTGeo::PHGeomIOTGeo() :
    _fGeom(NULL)
{
//  SplitLevel(0);
}

PHGeomIOTGeo::PHGeomIOTGeo(const PHGeomIOTGeo& geom) :
    _fGeom(NULL)
{
  if (geom.isValid())
    SetGeometry(geom._fGeom);
}

PHGeomIOTGeo::~PHGeomIOTGeo()
{
  Reset();
}

PHObject*
PHGeomIOTGeo::clone() const
{
  PHGeomIOTGeo * geo = new PHGeomIOTGeo(*this);
  return geo;
}

void
PHGeomIOTGeo::SetGeometry(const TGeoVolume * g)
{
  if (!g)
    {
      cout << __PRETTY_FUNCTION__ << " - Error - Invalid input" << endl;
      return;
    }

  if (_fGeom)
    delete _fGeom;

  _fGeom = static_cast< TGeoVolume *> (g->Clone());
}

TGeoVolume *
PHGeomIOTGeo::GetGeometryCopy() const
{
  if (not isValid()) return NULL;
  return static_cast< TGeoVolume *> (_fGeom->Clone());
}

/** identify Function from PHObject
 @param os Output Stream
 */
void
PHGeomIOTGeo::identify(std::ostream& os) const
{
  os << "PHGeomIOTGeo - ";
  if (_fGeom)
    os << " with geometry data " << _fGeom->GetName() << ": "
        << _fGeom->GetTitle();
  else
    os << "Empty";
  os << endl;
}

/// Clear Event
void
PHGeomIOTGeo::Reset()
{
  if (_fGeom)
    delete _fGeom;
  _fGeom = NULL;
}

/// isValid returns non zero if object contains vailid data
int
PHGeomIOTGeo::isValid() const
{
  if (_fGeom == NULL)
    return 0;
  if (_fGeom->IsZombie())
    return 0;
  return 1;
}
