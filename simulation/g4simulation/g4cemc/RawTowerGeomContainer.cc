#include "NewGeomContainer.h"
#include "NewGeom.h"

#include <cstdlib>
#include <iostream>

ClassImp(NewGeomContainer)

using namespace std;


NewGeomContainer::NewGeomContainer( RawTowerDefs::CalorimeterId caloid )
  {
    _caloid = caloid;
  }

NewGeomContainer::~NewGeomContainer()
{}

NewGeomContainer::ConstRange
NewGeomContainer::get_tower_geometries( void ) const
{
  return make_pair(_geoms.begin(), _geoms.end());
}


NewGeomContainer::Range
NewGeomContainer::get_tower_geometries( void )
{
  return make_pair(_geoms.begin(), _geoms.end());
}


NewGeomContainer::ConstIterator
NewGeomContainer::add_tower_geometry(NewGeom *geo)
{
  _geoms[geo->get_id()] = geo;
  return _geoms.find(geo->get_id());
}

NewGeom *
NewGeomContainer::get_tower_geometry(RawTowerDefs::keytype key)
{
  Iterator it = _geoms.find(key);
  if (it != _geoms.end())
    {
      return it->second;
    }
  return NULL;
}


int
NewGeomContainer::isValid() const
{
  return (!_geoms.empty());
}

void
NewGeomContainer::Reset()
{
  while (_geoms.begin() != _geoms.end())
    {
      delete _geoms.begin()->second;
      _geoms.erase(_geoms.begin());
    }
}

void
NewGeomContainer::identify(std::ostream& os) const
{
  os << "NewGeomContainer, number of tower geometries: " << size() << std::endl;
}
