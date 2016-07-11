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
#include <iostream>

using namespace std;

ClassImp(PHGeomTGeo);

PHGeomTGeo::PHGeomTGeo() :
    _fGeom(NULL)
{
  // TODO Auto-generated constructor stub

}

PHGeomTGeo::~PHGeomTGeo()
{
  // TODO Auto-generated destructor stub
}

void
PHGeomTGeo::SetGeometry(TGeoManager * g)
{
  if (!g)
    {
      cout << "dbFvtxAlignment::SetGeometry - Error - Invalid input" << endl;
      return;
    }

//  if (_fGeom)
//    {
//      cout << "dbFvtxAlignment::SetGeometry - Clean up the old geometry"
//          << endl;
//      delete _fGeom;
//    }

////  preserve the gGeoManager
//  TGeoManager * g_tmp = gGeoManager;
//  gGeoManager = NULL;
//
//  _fGeom = static_cast<TGeoManager *>(g->Clone("PHGeomTGeo_fGeom"));
//
//  // recover it
//  gGeoManager = g_tmp;

  _fGeom = g;

}

TGeoManager *
PHGeomTGeo::GetGeometry()
{
  return _fGeom;
}

/** identify Function from PHObject
 @param os Output Stream
 */
void
PHGeomTGeo::identify(std::ostream& os) const
{
  os << "PHGeomTGeo - ";
  if (_fGeom)
    os << " with geometry data " << _fGeom->GetName() << ": "
        << _fGeom->GetTitle();
  else
    os << "Empty";
  os << endl;
}

/// Clear Event
void
PHGeomTGeo::Reset()
{
  if (_fGeom)
    delete _fGeom;
  _fGeom = NULL;
}

/// isValid returns non zero if object contains vailid data
int
PHGeomTGeo::isValid() const
{
  if (_fGeom == NULL)
    return 0;
  if (_fGeom->IsZombie())
    return 0;
  return 1;
}
