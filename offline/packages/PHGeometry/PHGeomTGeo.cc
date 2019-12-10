// $Id: $

/*!
 * \file PHGeomTGeo.cc
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHGeomTGeo.h"

#include "TGeoManager.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

using namespace std;

PHGeomTGeo::PHGeomTGeo()
  : _fGeom(nullptr)
{
}

PHGeomTGeo::~PHGeomTGeo()
{
  ConsistencyCheck();
  if (_fGeom)
  {
    _fGeom->UnlockGeometry();
  }
    delete _fGeom;
}

void PHGeomTGeo::SetGeometry(TGeoManager* g)
{
  ConsistencyCheck();
  assert(_fGeom == nullptr);

  if (!g)
  {
    cout << __PRETTY_FUNCTION__ << " - Error - Invalid input" << endl;
    return;
  }

  _fGeom = g;
  _fGeom->LockGeometry();

  ConsistencyCheck();
}

TGeoManager*
PHGeomTGeo::GetGeometry()
{
  if (_fGeom == nullptr)
    return nullptr;

  ConsistencyCheck();

  if (_fGeom == gGeoManager)
    return _fGeom;
  else
  {
    return nullptr;
  }
}

/** identify Function from PHObject
 @param os Output Stream
 */
void PHGeomTGeo::identify(std::ostream& os) const
{
  os << "PHGeomTGeo - ";
  if (_fGeom)
    os << " with geometry data " << _fGeom->GetName() << ": "
       << _fGeom->GetTitle();
  else
    os << "Empty";
  os << endl;
  ConsistencyCheck();
}

/// Clear Event
void PHGeomTGeo::Reset()
{
  ConsistencyCheck();

  if (_fGeom)
  {
    _fGeom->UnlockGeometry();
    delete _fGeom;
  }
  _fGeom = nullptr;
}

/// isValid returns non zero if object contains vailid data
int PHGeomTGeo::isValid() const
{
  ConsistencyCheck();

  if (_fGeom == nullptr)
    return 0;
  if (_fGeom->IsZombie())
    return 0;
  return 1;
}

bool PHGeomTGeo::ConsistencyCheck() const
{
  if (_fGeom == nullptr)
    return true;  // uninitialized

  if (_fGeom == gGeoManager)
    return true;
  else
  {
    cout << __PRETTY_FUNCTION__
         << " - ERROR - gGeoManager is overridden by another TGeoManager. "
         << "Please avoid using multiple TGeoManager in processing. Stop the process."
         << endl;
    exit(1);
    return false;
  }
}
