#include "RawTowerGeomContainer_Cylinderv1.h"
#include "NewGeom.h"

#include <cstdlib>
#include <iostream>

ClassImp(RawTowerGeomContainer_Cylinderv1)

using namespace std;


RawTowerGeomContainer_Cylinderv1::RawTowerGeomContainer_Cylinderv1( RawTowerDefs::CalorimeterId caloid )
  {
    _caloid = caloid;
  }

RawTowerGeomContainer_Cylinderv1::~RawTowerGeomContainer_Cylinderv1()
{}

RawTowerGeomContainer_Cylinderv1::ConstRange
RawTowerGeomContainer_Cylinderv1::get_tower_geometries( void ) const
{
  return make_pair(_geoms.begin(), _geoms.end());
}


RawTowerGeomContainer_Cylinderv1::Range
RawTowerGeomContainer_Cylinderv1::get_tower_geometries( void )
{
  return make_pair(_geoms.begin(), _geoms.end());
}


RawTowerGeomContainer_Cylinderv1::ConstIterator
RawTowerGeomContainer_Cylinderv1::add_tower_geometry(NewGeom *geo)
{
  _geoms[geo->get_id()] = geo;
  return _geoms.find(geo->get_id());
}

NewGeom *
RawTowerGeomContainer_Cylinderv1::get_tower_geometry(RawTowerDefs::keytype key)
{
  Iterator it = _geoms.find(key);
  if (it != _geoms.end())
    {
      return it->second;
    }
  return NULL;
}


int
RawTowerGeomContainer_Cylinderv1::isValid() const
{
  return (!_geoms.empty());
}

void
RawTowerGeomContainer_Cylinderv1::Reset()
{
  while (_geoms.begin() != _geoms.end())
    {
      delete _geoms.begin()->second;
      _geoms.erase(_geoms.begin());
    }
}

void
RawTowerGeomContainer_Cylinderv1::identify(std::ostream& os) const
{
  os << "RawTowerGeomContainer_Cylinderv1, number of tower geometries: " << size() << std::endl;
}
