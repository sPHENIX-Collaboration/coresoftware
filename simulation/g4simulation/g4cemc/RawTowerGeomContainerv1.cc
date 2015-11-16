#include "RawTowerGeomContainerv1.h"
#include "RawTowerGeom.h"

#include <cstdlib>
#include <iostream>

ClassImp(RawTowerGeomContainerv1)

using namespace std;


RawTowerGeomContainerv1::RawTowerGeomContainerv1( RawTowerDefs::CalorimeterId caloid )
  {
    _caloid = caloid;
  }

RawTowerGeomContainerv1::~RawTowerGeomContainerv1()
{}

RawTowerGeomContainerv1::ConstRange
RawTowerGeomContainerv1::get_tower_geometries( void ) const
{
  return make_pair(_geoms.begin(), _geoms.end());
}


RawTowerGeomContainerv1::Range
RawTowerGeomContainerv1::get_tower_geometries( void )
{
  return make_pair(_geoms.begin(), _geoms.end());
}


RawTowerGeomContainerv1::ConstIterator
RawTowerGeomContainerv1::add_tower_geometry(RawTowerGeom *geo)
{
  _geoms[geo->get_id()] = geo;
  return _geoms.find(geo->get_id());
}

RawTowerGeom *
RawTowerGeomContainerv1::get_tower_geometry(RawTowerDefs::keytype key)
{
  Iterator it = _geoms.find(key);
  if (it != _geoms.end())
    {
      return it->second;
    }
  return NULL;
}


int
RawTowerGeomContainerv1::isValid() const
{
  return (!_geoms.empty());
}

void
RawTowerGeomContainerv1::Reset()
{
  while (_geoms.begin() != _geoms.end())
    {
      delete _geoms.begin()->second;
      _geoms.erase(_geoms.begin());
    }
}

void
RawTowerGeomContainerv1::identify(std::ostream& os) const
{
  os << "RawTowerGeomContainerv1, number of tower geometries: " << size() << std::endl;
}
