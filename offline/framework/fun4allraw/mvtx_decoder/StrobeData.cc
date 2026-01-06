#include "StrobeData.h"

#include <iostream>

///_________________________________________________________________
/// clear
void mvtx::StrobeData::clear()
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

