#include "EpdGeomV1.h"

#include <map>
#include <utility>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>


EpdGeomV1::EpdGeomV1()
{

}


EpdGeomV1::~EpdGeomV1() {};


unsigned int EpdGeomV1::get_arm_index(unsigned int key)
{
  unsigned int arm_index = key >> 20U;
  return arm_index;
}
unsigned int EpdGeomV1::get_r_index(unsigned int key)
{
  unsigned int arm_index = get_arm_index(key);
  unsigned int r_index = (key - (arm_index << 20U)) >> 10U;
  return r_index;
}
unsigned int EpdGeomV1::get_phi_index(unsigned int key)
{
  unsigned int arm_index = get_arm_index(key);
  unsigned int r_index = get_r_index(key);
  unsigned int phi_index = (key - (arm_index << 20U) - (r_index << 10U));
  return phi_index;
}



