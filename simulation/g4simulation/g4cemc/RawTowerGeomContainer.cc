#include "RawTowerGeomContainer.h"
#include "RawTowerGeom.h"

#include <cstdlib>
#include <iostream>

ClassImp(RawTowerGeomContainer)

using namespace std;


RawTowerGeomContainer::RawTowerGeomContainer( RawTowerDefs::CalorimeterId caloid )
  {
    _caloid = caloid;
  }

RawTowerGeomContainer::~RawTowerGeomContainer()
{}

RawTowerGeomContainer::ConstRange
RawTowerGeomContainer::get_tower_geometries( void ) const
{
  return make_pair(_geoms.begin(), _geoms.end());
}


RawTowerGeomContainer::Range
RawTowerGeomContainer::get_tower_geometries( void )
{
  return make_pair(_geoms.begin(), _geoms.end());
}


RawTowerGeomContainer::ConstIterator
RawTowerGeomContainer::add_tower_geometry(RawTowerGeom *geo)
{
  _geoms[geo->get_id()] = geo;
  return _geoms.find(geo->get_id());
}

RawTowerGeom *
RawTowerGeomContainer::get_tower_geometry(RawTowerDefs::keytype key)
{
  Iterator it = _geoms.find(key);
  if (it != _geoms.end())
    {
      return it->second;
    }
  return NULL;
}


int
RawTowerGeomContainer::isValid() const
{
  return (!_geoms.empty());
}

void
RawTowerGeomContainer::Reset()
{
  while (_geoms.begin() != _geoms.end())
    {
      delete _geoms.begin()->second;
      _geoms.erase(_geoms.begin());
    }
}

void
RawTowerGeomContainer::identify(std::ostream& os) const
{
  os << "RawTowerGeomContainer, number of tower geometries: " << size() << std::endl;
}
