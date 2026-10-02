#include "StrobeData.h"

///_________________________________________________________________
/// clear
void mvtx_offline::StrobeData::clear()
{
  ir.clear();
  hasCDW = false;
  calWord = {};

  for (auto&& hit : hit_vector)
  {
    delete hit;
  }
  hit_vector.clear();
}
