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
