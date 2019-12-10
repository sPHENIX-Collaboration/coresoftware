#include "PHG4HitDefs.h"

#include <tr1/functional>

namespace PHG4HitDefs
{

  int get_volume_id(const std::string & nodename)
  {
    return std::tr1::hash<std::string>()(nodename);
  }

}


