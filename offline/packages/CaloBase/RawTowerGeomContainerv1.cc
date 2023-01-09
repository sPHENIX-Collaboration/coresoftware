#include "RawTowerGeomContainerv1.h"

#include "RawTowerGeom.h"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

using namespace std;

RawTowerGeomContainerv1::RawTowerGeomContainerv1(RawTowerDefs::CalorimeterId caloid)
  : _caloid(caloid)
{
}

RawTowerGeomContainerv1::~RawTowerGeomContainerv1()
{
  Reset();  // make sure everything is deleted
}

RawTowerGeomContainerv1::ConstRange
RawTowerGeomContainerv1::get_tower_geometries() const
{
  return make_pair<ConstIterator, ConstIterator>(_geoms.begin(), _geoms.end());
}

RawTowerGeomContainerv1::Range
RawTowerGeomContainerv1::get_tower_geometries()
{
  return make_pair<Iterator, Iterator>(_geoms.begin(), _geoms.end());
}

RawTowerGeomContainerv1::ConstIterator
RawTowerGeomContainerv1::add_tower_geometry(RawTowerGeom* geo)
{
  assert(geo);

  if (RawTowerDefs::decode_caloid(geo->get_id()) != get_calorimeter_id())
  {
    cout << "RawTowerGeomContainerv1::add_tower_geometry - Fatal Error - "
            "attempting to add tower geometry with id = "
         << geo->get_id()
         << " with CaloID = " << RawTowerDefs::decode_caloid(geo->get_id())
         << " to this container of CaloID = " << get_calorimeter_id() << ".";
    geo->identify(cout);
    exit(2);
  }

  Iterator it = _geoms.find(geo->get_id());
  if (it != _geoms.end())
  {
    cout
        << "RawTowerGeomContainerv1::add_tower_geometry - WARNING - replace tower geometry for tower #"
        << geo->get_id() << ". This Old tower will be deleted: ";
    it->second->identify(cout);

    delete it->second;
    _geoms.erase(it);
  }

  _geoms[geo->get_id()] = geo;
  return _geoms.find(geo->get_id());
}

RawTowerGeom*
RawTowerGeomContainerv1::get_tower_geometry(RawTowerDefs::keytype key)
{
  Iterator it = _geoms.find(key);
  if (it != _geoms.end())
  {
    return it->second;
  }
  return nullptr;
}

int RawTowerGeomContainerv1::isValid() const
{
  return (!_geoms.empty());
}

void RawTowerGeomContainerv1::Reset()
{
  while (_geoms.begin() != _geoms.end())
  {
    delete _geoms.begin()->second;
    _geoms.erase(_geoms.begin());
  }
}

void RawTowerGeomContainerv1::identify(std::ostream& os) const
{
  os << "RawTowerGeomContainerv1, number of tower geometries: " << size()
     << std::endl;
}
