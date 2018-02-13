#include "RawTowerGeomContainer_Planev1.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

RawTowerGeomContainer_Planev1::RawTowerGeomContainer_Planev1(RawTowerDefs::CalorimeterId caloid) :
  RawTowerGeomContainerv1(caloid),
  _center_x(NAN),
  _center_y(NAN),
  _center_z(NAN),
  _theta(NAN),
  _phi(NAN)
{
  return;
}


void
RawTowerGeomContainer_Planev1::Reset()
{

  RawTowerGeomContainerv1::Reset();
}


