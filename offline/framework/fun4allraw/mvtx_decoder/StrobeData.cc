#include "StrobeData.h"
#include <iostream>

using namespace mvtx;

///_________________________________________________________________
/// Destructor
StrobeData::~StrobeData()
{
}

///_________________________________________________________________
/// clear
void StrobeData::clear()
{
  ir.clear();
  hasCDW = false;
  calWord = {};

  for ( auto&& hit : hit_vector )
  {
    delete hit;
  }
  hit_vector.clear();
}

