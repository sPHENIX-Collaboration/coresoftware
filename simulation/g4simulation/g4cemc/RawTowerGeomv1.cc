#include "NewGeomv1.h"

#include <g4main/PHG4Utils.h>

#include <iostream>
#include <algorithm>

#include <cmath>
#include <map>

using namespace std;

ClassImp(NewGeomv1)

NewGeomv1::NewGeomv1(RawTowerDefs::keytype id) :
  _towerid(id),
  _center_x(0),
  _center_y(0),
  _center_z(0),
  _size_x(0),
  _size_y(0),
  _size_z(0)
{}

NewGeomv1::~NewGeomv1()
{}


double NewGeomv1::get_center_radius() const
{
  return sqrt( _center_x * _center_x +
	       _center_y * _center_y +
	       _center_z * _center_z );
}

double NewGeomv1::get_eta() const
{
  std::pair<double,double> etaphi;
  etaphi = PHG4Utils::get_etaphi( _center_x , _center_y , _center_z );

  return etaphi.first;
}

double NewGeomv1::get_phi() const
{
  std::pair<double,double> etaphi;
  etaphi = PHG4Utils::get_etaphi( _center_x , _center_y , _center_z );

  return etaphi.second;
}


void NewGeomv1::identify(std::ostream& os) const
{
  std::cout << "NewGeomv1:  x: " << get_center_x() << "  y: " << get_center_y() << "  z: " << get_center_z()
	    << "\n           dx: " << get_size_x() << " dy: " << get_size_y() << " dz: " << get_size_z() << std::endl;
}
